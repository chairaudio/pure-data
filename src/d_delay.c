/* Copyright (c) 1997-1999 Miller Puckette.
* For information on usage and redistribution, and for a DISCLAIMER OF ALL
* WARRANTIES, see the file, "LICENSE.txt," in this distribution.  */

/*  send~, delread~, throw~, catch~ */

#include "m_pd.h"
#include <string.h>
extern int ugen_getsortno(void);

#define DEFDELVS 64             /* LATER get this from canvas at DSP time */
static const int delread_zero = 0;    /* four bytes of zero for delread~, vd~/delread4~*/

/* ----------------------------- delwrite~ ----------------------------- */
static t_class *sigdelwrite_class;

typedef struct delwritectl
{
    int c_n;
    t_sample *c_vec;
    int c_phase;
} t_delwritectl;

typedef struct _sigdelwrite
{
    t_object x_obj;
    t_symbol *x_sym;
    t_float x_deltime;  /* delay in msec (added by Mathieu Bouchard) */
    t_delwritectl x_cspace;
    int x_sortno;   /* DSP sort number at which this was last put on chain */
    int x_rsortno;  /* DSP sort # for first delread or write in chain */
    int x_vecsize;  /* vector size for delread~ to use */
    t_float x_f;
} t_sigdelwrite;

#define XTRASAMPS 4
#define SAMPBLK 4

static void sigdelwrite_updatesr(t_sigdelwrite *x, t_float sr) /* added by Mathieu Bouchard */
{
    int nsamps = x->x_deltime * sr * (t_float)(0.001f);
    if (nsamps < 1) nsamps = 1;
    nsamps += ((- nsamps) & (SAMPBLK - 1));
    nsamps += DEFDELVS;
    if (x->x_cspace.c_n != nsamps)
    {
        x->x_cspace.c_vec = (t_sample *)resizebytes(x->x_cspace.c_vec,
            (x->x_cspace.c_n + XTRASAMPS) * sizeof(t_sample),
            (nsamps + XTRASAMPS) * sizeof(t_sample));
        x->x_cspace.c_n = nsamps;
        x->x_cspace.c_phase = XTRASAMPS;
    }
}

static void sigdelwrite_clear (t_sigdelwrite *x) /* added by Orm Finnendahl */
{
  if (x->x_cspace.c_n > 0)
    memset(x->x_cspace.c_vec, 0, sizeof(t_sample)*(x->x_cspace.c_n + XTRASAMPS));
}


    /* routine to check that all delwrites/delreads/vds have same vecsize */
static void sigdelwrite_checkvecsize(t_sigdelwrite *x, int vecsize)
{
    if (x->x_rsortno != ugen_getsortno())
    {
        x->x_vecsize = vecsize;
        x->x_rsortno = ugen_getsortno();
    }
    /*
        LATER this should really check sample rate and blocking, once that is
        supported.  Probably we don't actually care about vecsize.
        For now just suppress this check. */
#if 0
    else if (vecsize != x->x_vecsize)
        pd_error(x, "delread/delwrite/vd vector size mismatch");
#endif
}

static void *sigdelwrite_new(t_symbol *s, t_floatarg msec)
{
    t_sigdelwrite *x = (t_sigdelwrite *)pd_new(sigdelwrite_class);
    if (!*s->s_name) s = gensym("delwrite~");
    pd_bind(&x->x_obj.ob_pd, s);
    x->x_sym = s;
    x->x_deltime = msec;
    x->x_cspace.c_n = 0;
    x->x_cspace.c_vec = getbytes(XTRASAMPS * sizeof(t_sample));
    x->x_sortno = 0;
    x->x_vecsize = 0;
    x->x_f = 0;
    return (x);
}

static t_int *sigdelwrite_perform(t_int *w)
{
    t_sample *in = (t_sample *)(w[1]);
    t_delwritectl *c = (t_delwritectl *)(w[2]);
    int n = (int)(w[3]);
    int phase = c->c_phase, nsamps = c->c_n;
    t_sample *vp = c->c_vec, *bp = vp + phase, *ep = vp + (c->c_n + XTRASAMPS);
    phase += n;

    while (n--)
    {
        t_sample f = *in++;
        if (PD_BIGORSMALL(f))
            f = 0;
        *bp++ = f;
        if (bp == ep)
        {
            vp[0] = ep[-4];
            vp[1] = ep[-3];
            vp[2] = ep[-2];
            vp[3] = ep[-1];
            bp = vp + XTRASAMPS;
            phase -= nsamps;
        }
    }
    c->c_phase = phase;
    return (w+4);
}

static void sigdelwrite_dsp(t_sigdelwrite *x, t_signal **sp)
{
    dsp_add(sigdelwrite_perform, 3, sp[0]->s_vec, &x->x_cspace, (t_int)sp[0]->s_n);
    x->x_sortno = ugen_getsortno();
    sigdelwrite_checkvecsize(x, sp[0]->s_n);
    sigdelwrite_updatesr(x, sp[0]->s_sr);
}

static void sigdelwrite_free(t_sigdelwrite *x)
{
    pd_unbind(&x->x_obj.ob_pd, x->x_sym);
    freebytes(x->x_cspace.c_vec,
        (x->x_cspace.c_n + XTRASAMPS) * sizeof(t_sample));
}

static void sigdelwrite_setup(void)
{
    sigdelwrite_class = class_new(gensym("delwrite~"),
        (t_newmethod)sigdelwrite_new, (t_method)sigdelwrite_free,
        sizeof(t_sigdelwrite), 0, A_DEFSYM, A_DEFFLOAT, 0);
    CLASS_MAINSIGNALIN(sigdelwrite_class, t_sigdelwrite, x_f);
    class_addmethod(sigdelwrite_class, (t_method)sigdelwrite_dsp,
        gensym("dsp"), A_CANT, 0);
    class_addmethod(sigdelwrite_class, (t_method)sigdelwrite_clear,
                    gensym("clear"), 0);
}

/* ----------------------------- delread~ ----------------------------- */
static t_class *sigdelread_class;

typedef struct _sigdelread
{
    t_object x_obj;
    t_symbol *x_sym;
    t_float x_deltime;  /* delay in msec */
    int x_delsamps;     /* delay in samples */
    t_float x_sr;       /* samples per msec */
    t_float x_n;        /* vector size */
    int x_zerodel;      /* 0 or vecsize depending on read/write order */
} t_sigdelread;

static void sigdelread_float(t_sigdelread *x, t_float f);

static void *sigdelread_new(t_symbol *s, t_floatarg f)
{
    t_sigdelread *x = (t_sigdelread *)pd_new(sigdelread_class);
    x->x_sym = s;
    x->x_sr = 1;
    x->x_n = 1;
    x->x_zerodel = 0;
    sigdelread_float(x, f);
    outlet_new(&x->x_obj, &s_signal);
    return (x);
}

static void sigdelread_float(t_sigdelread *x, t_float f)
{
    t_sigdelwrite *delwriter =
        (t_sigdelwrite *)pd_findbyclass(x->x_sym, sigdelwrite_class);
    x->x_deltime = f;
    if (delwriter)
    {
        x->x_delsamps = (int)(0.5 + x->x_sr * x->x_deltime)
            + x->x_n - x->x_zerodel;
        if (x->x_delsamps < x->x_n) x->x_delsamps = x->x_n;
        else if (x->x_delsamps > delwriter->x_cspace.c_n)
            x->x_delsamps = delwriter->x_cspace.c_n;
    }
}

static t_int *sigdelread_perform(t_int *w)
{
    t_sample *out = (t_sample *)(w[1]);
    t_delwritectl *c = (t_delwritectl *)(w[2]);
    int delsamps = *(int *)(w[3]);
    int n = (int)(w[4]);
    int phase = c->c_phase - delsamps, nsamps = c->c_n;
    t_sample *vp = c->c_vec, *bp, *ep = vp + (c->c_n + XTRASAMPS);
    if (phase < 0) phase += nsamps;
    bp = vp + phase;

    while (n--)
    {
        *out++ = *bp++;
        if (bp == ep) bp -= nsamps;
    }
    return (w+5);
}

static void sigdelread_dsp(t_sigdelread *x, t_signal **sp)
{
    t_sigdelwrite *delwriter =
        (t_sigdelwrite *)pd_findbyclass(x->x_sym, sigdelwrite_class);
    x->x_sr = sp[0]->s_sr * 0.001;
    x->x_n = sp[0]->s_n;
    if (delwriter)
    {
        sigdelwrite_updatesr(delwriter, sp[0]->s_sr);
        sigdelwrite_checkvecsize(delwriter, sp[0]->s_n);
        x->x_zerodel = (delwriter->x_sortno == ugen_getsortno() ?
            0 : delwriter->x_vecsize);
        sigdelread_float(x, x->x_deltime);
        dsp_add(sigdelread_perform, 4,
            sp[0]->s_vec, &delwriter->x_cspace, &x->x_delsamps, (t_int)sp[0]->s_n);
        /* check block size - but only if delwriter has been initialized */
        if (delwriter->x_cspace.c_n > 0 && sp[0]->s_n > delwriter->x_cspace.c_n)
            pd_error(x, "delread~ %s: blocksize larger than delwrite~ buffer", x->x_sym->s_name);
    }
    else if (*x->x_sym->s_name)
        pd_error(x, "delread~: %s: no such delwrite~",x->x_sym->s_name);
}

static void sigdelread_setup(void)
{
    sigdelread_class = class_new(gensym("delread~"),
        (t_newmethod)sigdelread_new, 0,
        sizeof(t_sigdelread), 0, A_DEFSYM, A_DEFFLOAT, 0);
    class_addmethod(sigdelread_class, (t_method)sigdelread_dsp,
        gensym("dsp"), A_CANT, 0);
    class_addfloat(sigdelread_class, (t_method)sigdelread_float);
}


/* ----------------------------- vd~ / delread4~ ----------------------------- */
static t_class *sigvd_class;

typedef struct _sigvd
{
    t_object x_obj;
    t_symbol *x_sym;
    t_float x_sr;       /* samples per msec */
    int x_zerodel;      /* 0 or vecsize depending on read/write order */
    t_float x_f;
} t_sigvd;

static void *sigvd_new(t_symbol *s)
{
    t_sigvd *x = (t_sigvd *)pd_new(sigvd_class);
    x->x_sym = s;
    x->x_sr = 1;
    x->x_zerodel = 0;
    outlet_new(&x->x_obj, &s_signal);
    x->x_f = 0;
    return (x);
}

static t_int *sigvd_perform(t_int *w)
{
    t_sample *in = (t_sample *)(w[1]);
    t_sample *out = (t_sample *)(w[2]);
    t_delwritectl *ctl = (t_delwritectl *)(w[3]);
    t_sigvd *x = (t_sigvd *)(w[4]);
    int n = (int)(w[5]);

    int nsamps = ctl->c_n;
    t_sample limit = nsamps - n;
    t_sample fn = n-1;
    t_sample *vp = ctl->c_vec, *bp, *wp = vp + ctl->c_phase;
    t_sample zerodel = x->x_zerodel;
    if (limit < 0) /* blocksize is larger than delread~ buffer size */
    {
        while (n--)
            *out++ = 0;
        return (w+6);
    }
    while (n--)
    {
        t_sample delsamps = x->x_sr * *in++ - zerodel, frac;
        int idelsamps;
        t_sample a, b, c, d, cminusb;
        if (!(delsamps >= 1.00001f))    /* too small or NAN */
            delsamps = 1.00001f;
        if (delsamps > limit)           /* too big */
            delsamps = limit;
        delsamps += fn;
        fn = fn - 1.0f;
        idelsamps = delsamps;
        frac = delsamps - (t_sample)idelsamps;
        bp = wp - idelsamps;
        if (bp < vp + XTRASAMPS) bp += nsamps;
        d = bp[-3];
        c = bp[-2];
        b = bp[-1];
        a = bp[0];
        cminusb = c-b;
        *out++ = b + frac * (
            cminusb - 0.1666667f * (1.-frac) * (
                (d - a - 3.0f * cminusb) * frac + (d + 2.0f*a - 3.0f*b)
            )
        );
    }
    return (w+6);
}

static void sigvd_dsp(t_sigvd *x, t_signal **sp)
{
    t_sigdelwrite *delwriter =
        (t_sigdelwrite *)pd_findbyclass(x->x_sym, sigdelwrite_class);
    x->x_sr = sp[0]->s_sr * 0.001;
    if (delwriter)
    {
        sigdelwrite_checkvecsize(delwriter, sp[0]->s_n);
        x->x_zerodel = (delwriter->x_sortno == ugen_getsortno() ?
            0 : delwriter->x_vecsize);
        dsp_add(sigvd_perform, 5,
            sp[0]->s_vec, sp[1]->s_vec,
                &delwriter->x_cspace, x, (t_int)sp[0]->s_n);
        /* check block size - but only if delwriter has been initialized */
        if (delwriter->x_cspace.c_n > 0 && sp[0]->s_n > delwriter->x_cspace.c_n)
            pd_error(x, "delread4~ %s: blocksize larger than delwrite~ buffer", x->x_sym->s_name);
    }
    else if (*x->x_sym->s_name)
        pd_error(x, "delread4~: %s: no such delwrite~",x->x_sym->s_name);
}

static void sigvd_setup(void)
{
    sigvd_class = class_new(gensym("delread4~"), (t_newmethod)sigvd_new, 0,
        sizeof(t_sigvd), 0, A_DEFSYM, 0);
    class_addcreator((t_newmethod)sigvd_new, gensym("vd~"), A_DEFSYM, 0);
    class_addmethod(sigvd_class, (t_method)sigvd_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(sigvd_class, t_sigvd, x_f);
}

/* ----------------------------- delreadsinc~ ----------------------------- */
// the following defines the quality of interpolation
// the values have been chosen to provide high quality with 48Khz sample rate and 20kHz pass band
// stop band rejection is 23dB at 24kHz and 87dB at 28kHz wich would be the first frequency
// to back-alias into the 20kHz range
//
#define STEPS_ZC 4096 // table values per zero crossing this is accurate for 24 bits precision
#define HALF_N_ZC 11 // zero crossings on pos and neg axis 
#define LP_SCALE 0.85  // defines cuttof frequency of sinc lowpass in normalized frequency
#define N_SMP 2*HALF_N_ZC+1 // number of samples for convolution
#define SINC_LEN (STEPS_ZC * (HALF_N_ZC+1)) // sin function table length
#include <math.h> // for sinc table generation

static t_class *sigvdsinc_class;

typedef struct _sigvdsinc
{
    t_object x_obj;
    t_symbol *x_sym;
    t_float x_sr;       /* samples per msec */
    int x_zerodel;      /* 0 or vecsize depending on read/write order */
    t_float x_f;
      // sinc table is held in member variable
    t_sample sinc_array[SINC_LEN];
    // derivative of sinc funtion for interpolation of sinc function table
    t_sample sinc_diff_array[SINC_LEN];
} t_sigvdsinc;

void sigvdsinc_initialize_sinc_table(t_sample* sinc_array, t_sample* sinc_diff_array) {
  
  sinc_array[0] = LP_SCALE;
  sinc_diff_array[SINC_LEN] = 0;
  
  float a0 = 0.35875;
  float a1 = 0.48829;
  float a2 = 0.14128;
  float a3 = 0.01168;

  for (int i=1; i<SINC_LEN; i++) {
    float idx = 0.5*(M_PI*i/SINC_LEN - M_PI);
    // four term blackmanharris after https://www.mathworks.com/help/signal/ref/blackmanharris.html
    float blackmanharris= a0 - a1*cos(2*idx) + a2*cos(4*idx) - a3*cos(6*idx);
    // blackmanharris windowed sinc function
    sinc_array[i] = sin(LP_SCALE*M_PI*(float)i/STEPS_ZC)/(LP_SCALE*M_PI*(float)i/STEPS_ZC) * blackmanharris * LP_SCALE;
    // calculate derivative for table interpolation
    sinc_diff_array[i-1] = sinc_array[i]-sinc_array[i-1];
  }
}

t_sample  sigvdsinc_interpolate(t_sample* samples, t_sample fraction, t_sample* sinc_array, t_sample* sinc_diff_array){
  // samples must be of odd length and should be at least of size 9
  t_sample y = 0;
  t_sample exact_idx = 0;
  int sincf_idx = 0;
  t_sample rest_idx = 0;

  // if fraction is negative, take anti-fraction and shift all samples to the right
  int s = 0;
  if(fraction<0){
    fraction=fraction+1;
    s = 1;
  }

  exact_idx = (fraction)*(t_sample)STEPS_ZC;
  sincf_idx = roundf(exact_idx);
  rest_idx = exact_idx-sincf_idx;
    // midpoint
  y = y+samples[HALF_N_ZC-s]*(sinc_array[sincf_idx]+rest_idx*sinc_diff_array[sincf_idx]);

  for(int i=1; i<=HALF_N_ZC; i++) {
    // negative half
    exact_idx = ((t_sample)i+fraction)*(t_sample)STEPS_ZC;
    sincf_idx = roundf(exact_idx);
    rest_idx = exact_idx-sincf_idx;
    y = y+samples[HALF_N_ZC-i-s]*(sinc_array[sincf_idx]+rest_idx*sinc_diff_array[sincf_idx]);
    // positive half
    exact_idx = ((t_sample)i-fraction)*(t_sample)STEPS_ZC;
    sincf_idx = roundf(exact_idx);
    rest_idx = exact_idx-sincf_idx;
    y = y+samples[HALF_N_ZC+i-s]*(sinc_array[sincf_idx]+rest_idx*sinc_diff_array[sincf_idx]);
  }
  return y;
}

static void *sigvdsinc_new(t_symbol *s)
{
    t_sigvdsinc *x = (t_sigvdsinc *)pd_new(sigvdsinc_class);
    x->x_sym = s;
    x->x_sr = 1;
    x->x_zerodel = 0;
    outlet_new(&x->x_obj, &s_signal);
    x->x_f = 0;
    sigvdsinc_initialize_sinc_table(x->sinc_array, x->sinc_diff_array);
    return (x);
}

static t_int *sigvdsinc_perform(t_int *w)
{
    t_sample *in = (t_sample *)(w[1]);
    t_sample *out = (t_sample *)(w[2]);
    t_delwritectl *ctl = (t_delwritectl *)(w[3]);
    t_sigvdsinc *x = (t_sigvdsinc *)(w[4]);
    int n = (int)(w[5]);

    int nsamps = ctl->c_n;
    t_sample limit = nsamps - (n+HALF_N_ZC+1);
    t_sample fn = n-1;
    t_sample *vp = ctl->c_vec;
    t_sample *bp;
    t_sample *wp = vp + ctl->c_phase;
    t_sample zerodel = x->x_zerodel;

    if (limit < 0 ) /* blocksize is larger than delread~ buffer size */
    {
        while (n--)
            *out++ = 0;
        return (w+6);
    }
    t_sample samples[N_SMP];
    while (n--)
    {
        t_sample frac;
        t_sample delsamps = x->x_sr * *in++ - zerodel;
        delsamps = delsamps +1;
        int idelsamps;
        if (delsamps < HALF_N_ZC)   /* too small or NAN */
            delsamps = HALF_N_ZC;
        if (delsamps > limit)           /* too big */
            delsamps = limit;
        delsamps += fn;
        fn = fn - 1.0f;
        idelsamps = roundf(delsamps);
        bp = wp - idelsamps;
        if (bp < vp + XTRASAMPS) bp += nsamps;

        frac = delsamps - (t_sample)idelsamps;

        if(frac == 0){
            *out++ = *bp;
        }        
        else{
            samples[HALF_N_ZC] = *bp;

            for(int j=1; j<=HALF_N_ZC; j++)
            {
              t_sample* right = bp+j;
              if(right< vp + XTRASAMPS) {
                right += nsamps;
              }
              else if(right > vp + nsamps + XTRASAMPS){
                right -=nsamps;
              }
              samples[HALF_N_ZC+j] = *right;

              t_sample* left = bp-j;
              if(left< vp + XTRASAMPS) left += nsamps;
              samples[HALF_N_ZC-j] = *left;
            }

            *out++ = sigvdsinc_interpolate( samples, -frac, x->sinc_array, x->sinc_diff_array);
        }
    }
    return (w+6);
}

static void sigvdsinc_dsp(t_sigvdsinc *x, t_signal **sp)
{
    t_sigdelwrite *delwriter =
        (t_sigdelwrite *)pd_findbyclass(x->x_sym, sigdelwrite_class);
    x->x_sr = sp[0]->s_sr * 0.001;
    if (delwriter)
    {
        sigdelwrite_checkvecsize(delwriter, sp[0]->s_n);
        x->x_zerodel = (delwriter->x_sortno == ugen_getsortno() ?
            0 : delwriter->x_vecsize);
        dsp_add(sigvdsinc_perform, 5,
            sp[0]->s_vec, sp[1]->s_vec,
                &delwriter->x_cspace, x, (t_int)sp[0]->s_n);
        /* check block size - but only if delwriter has been initialized */
        if (delwriter->x_cspace.c_n > 0 && sp[0]->s_n > delwriter->x_cspace.c_n)
            pd_error(x, "delreadsinc~ %s: blocksize larger than delwrite~ buffer", x->x_sym->s_name);
    }
    else if (*x->x_sym->s_name)
        pd_error(x, "delreadsinc~: %s: no such delwrite~",x->x_sym->s_name);
}

static void sigvdsinc_setup(void)
{
    sigvdsinc_class = class_new(gensym("delreadsinc~"), (t_newmethod)sigvdsinc_new, 0,
        sizeof(t_sigvdsinc), 0, A_DEFSYM, 0);
    class_addmethod(sigvdsinc_class, (t_method)sigvdsinc_dsp, gensym("dsp"), A_CANT, 0);
    CLASS_MAINSIGNALIN(sigvdsinc_class, t_sigvdsinc, x_f);
}

/* ----------------------- global setup routine ---------------- */

void d_delay_setup(void)
{
    sigdelwrite_setup();
    sigdelread_setup();
    sigvd_setup();
    sigvdsinc_setup();
}

