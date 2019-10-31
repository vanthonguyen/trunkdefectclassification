#ifndef BSPLINES_H
#define BSPLINES_H

#include <iostream>
#include <fstream>

#include <math.h>
#include <ctime>

//#include <gsl/gsl_bspline.h> 
#include <gsl/gsl_errno.h>

#include <gsl/gsl_spline.h>
#include <gsl/gsl_bspline.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_statistics.h>

#include "DGtal/base/Common.h"
#include "DGtal/helpers/StdDefs.h"


//#include "TubeAnalyseHelper.h"

using namespace DGtal;


class BSplines{
public:
        static std::vector<Z3i::RealPoint> bsplines(std::vector<Z3i::RealPoint> points, int nbCoeffs, const std::vector<double> &weights){
        std::vector<Z3i::RealPoint> lines;
        /* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
        int n = points.size();
        size_t nbBreak = nbCoeffs - 2;
        gsl_bspline_workspace *bw;
        gsl_vector *BX, *BY, *BZ;
        gsl_matrix *covX, *covY, *covZ;
        gsl_vector *cX, *cY, *cZ;
        gsl_vector *xs, *ys, *zs, *ts;
        gsl_vector *wx, *wy, *wz;
        //gsl_multifit_linear_workspace *mw;
       
        //breakpoint
        //gsl_vector *bp;

        xs = gsl_vector_alloc(n);
        ys = gsl_vector_alloc(n);
        zs = gsl_vector_alloc(n);
        ts = gsl_vector_alloc(n);

        wx = gsl_vector_alloc(n);
        wy = gsl_vector_alloc(n);
        wz = gsl_vector_alloc(n);


        cX = gsl_vector_alloc(nbCoeffs);
        cY = gsl_vector_alloc(nbCoeffs);
        cZ = gsl_vector_alloc(nbCoeffs);



        //bp = gsl_vector_alloc(nbBreak);

        covX = gsl_matrix_alloc(nbCoeffs, nbCoeffs);
        covY = gsl_matrix_alloc(nbCoeffs, nbCoeffs);
        covZ = gsl_matrix_alloc(nbCoeffs, nbCoeffs);

        //cubic
        bw = gsl_bspline_alloc(4, nbBreak);
        BX = gsl_vector_alloc(nbCoeffs);
        BY = gsl_vector_alloc(nbCoeffs);
        BZ = gsl_vector_alloc(nbCoeffs);

        for(unsigned int i = 0; i< points.size(); i++){
            gsl_vector_set(xs, i, points[i][0]);
            gsl_vector_set(ys, i, points[i][1]);
            gsl_vector_set(zs, i, points[i][2]);

            gsl_vector_set(wx, i, weights[i]);
            gsl_vector_set(wy, i, weights[i]);
            gsl_vector_set(wz, i, weights[i]);

            gsl_vector_set(ts, i, i);
//std::cout<< i << "  "<<points[i][2]<<std::endl;
        }
        
       // for(unsigned int i = 0; i< breakPoints.size(); i++){
       //     gsl_vector_set(bp, i, (double)breakPoints[i]);
       // }
        gsl_bspline_knots_uniform(0.0, n, bw);
        //gsl_bspline_knots(bp, bw);
        bsplines(xs, ts, wx,  bw, BX, cX, covX, nbCoeffs, n);
        bsplines(ys, ts, wy, bw, BY, cY, covY, nbCoeffs, n);
        bsplines(zs, ts, wz, bw, BZ, cZ, covZ, nbCoeffs, n);
        double xi, yi, zi, err;
        double dist = 0.05;
        bool first = true;
        Z3i::RealPoint lastAddedPoint;
        /**
         */
        for (double ti = 0.1; ti < n; ti += 0.1){
            gsl_bspline_eval(ti, BX, bw);
            gsl_multifit_linear_est(BX, cX, covX, &xi, &err);

            gsl_bspline_eval(ti, BY, bw);
            gsl_multifit_linear_est(BY, cY, covY, &yi, &err);

            gsl_bspline_eval(ti, BZ, bw);
            gsl_multifit_linear_est(BZ, cZ, covZ, &zi, &err);
            
            Z3i::RealPoint p(xi, yi, zi);
            if(first ||  (p - lastAddedPoint).norm() >= dist ){
                lines.push_back(p);
                lastAddedPoint = p;
                first = false;
            }
        }

        //free
        gsl_bspline_free(bw);

        gsl_vector_free(ts);
        gsl_vector_free(xs);
        gsl_vector_free(ys);
        gsl_vector_free(zs);

        gsl_vector_free(wx);
        gsl_vector_free(wy);
        gsl_vector_free(wz);

        gsl_vector_free(cX);
        gsl_vector_free(cY);
        gsl_vector_free(cZ);

        gsl_vector_free(BX);
        gsl_vector_free(BY);
        gsl_vector_free(BZ);

        //gsl_vector_free(bp);

        gsl_matrix_free(covX);
        gsl_matrix_free(covY);
        gsl_matrix_free(covZ);

        return lines;

    }

    static std::vector<Z3i::RealPoint> bsplines(std::vector<Z3i::RealPoint> points, int nbCoeffs){
        std::vector<Z3i::RealPoint> lines;
        /* nbreak = ncoeffs + 2 - k = ncoeffs - 2 since k = 4 */
        int n = points.size();
        size_t nbBreak = nbCoeffs - 2;
        gsl_bspline_workspace *bw;
        gsl_vector *BX, *BY, *BZ;
        gsl_matrix *covX, *covY, *covZ;
        gsl_vector *cX, *cY, *cZ;
        gsl_vector *xs, *ys, *zs, *ts;
        gsl_vector *wx, *wy, *wz;
        //gsl_multifit_linear_workspace *mw;
       
        //breakpoint
        //gsl_vector *bp;

        xs = gsl_vector_alloc(n);
        ys = gsl_vector_alloc(n);
        zs = gsl_vector_alloc(n);
        ts = gsl_vector_alloc(n);

        wx = gsl_vector_alloc(n);
        wy = gsl_vector_alloc(n);
        wz = gsl_vector_alloc(n);


        cX = gsl_vector_alloc(nbCoeffs);
        cY = gsl_vector_alloc(nbCoeffs);
        cZ = gsl_vector_alloc(nbCoeffs);



        //bp = gsl_vector_alloc(nbBreak);

        covX = gsl_matrix_alloc(nbCoeffs, nbCoeffs);
        covY = gsl_matrix_alloc(nbCoeffs, nbCoeffs);
        covZ = gsl_matrix_alloc(nbCoeffs, nbCoeffs);

        bw = gsl_bspline_alloc(4, nbBreak);
        BX = gsl_vector_alloc(nbCoeffs);
        BY = gsl_vector_alloc(nbCoeffs);
        BZ = gsl_vector_alloc(nbCoeffs);

        for(unsigned int i = 0; i< points.size(); i++){
            gsl_vector_set(xs, i, points[i][0]);
            gsl_vector_set(ys, i, points[i][1]);
            gsl_vector_set(zs, i, points[i][2]);

            //gsl_vector_set(wx, i, 1 / points[i][0] / points[i][0]);
            //gsl_vector_set(wy, i, 1 / points[i][1] / points[i][1]);
            //gsl_vector_set(wz, i, 1 / points[i][2] / points[i][2]);

            gsl_vector_set(wx, i, 1);
            gsl_vector_set(wy, i, 1);
            gsl_vector_set(wz, i, 1);
            gsl_vector_set(ts, i, i);
//std::cout<< i << "  "<<points[i][2]<<std::endl;
        }
        
       // for(unsigned int i = 0; i< breakPoints.size(); i++){
       //     gsl_vector_set(bp, i, (double)breakPoints[i]);
       // }
        gsl_bspline_knots_uniform(0.0, n, bw);
        //gsl_bspline_knots(bp, bw);
        bsplines(xs, ts, wx,  bw, BX, cX, covX, nbCoeffs, n);
        bsplines(ys, ts, wy, bw, BY, cY, covY, nbCoeffs, n);
        bsplines(zs, ts, wz, bw, BZ, cZ, covZ, nbCoeffs, n);
        double xi, yi, zi, err;
        double dist = 0.05;
        bool first = true;
        Z3i::RealPoint lastAddedPoint;
        for (double ti = 0.1; ti < n; ti += 0.1){
            gsl_bspline_eval(ti, BX, bw);
            gsl_multifit_linear_est(BX, cX, covX, &xi, &err);

            gsl_bspline_eval(ti, BY, bw);
            gsl_multifit_linear_est(BY, cY, covY, &yi, &err);

            gsl_bspline_eval(ti, BZ, bw);
            gsl_multifit_linear_est(BZ, cZ, covZ, &zi, &err);
            
            Z3i::RealPoint p(xi, yi, zi);
            if(first ||  (p - lastAddedPoint).norm() >= dist ){
                lines.push_back(p);
                lastAddedPoint = p;
                first = false;
            }
        }

        //free
        gsl_bspline_free(bw);

        gsl_vector_free(ts);
        gsl_vector_free(xs);
        gsl_vector_free(ys);
        gsl_vector_free(zs);

        gsl_vector_free(wx);
        gsl_vector_free(wy);
        gsl_vector_free(wz);

        gsl_vector_free(cX);
        gsl_vector_free(cY);
        gsl_vector_free(cZ);

        gsl_vector_free(BX);
        gsl_vector_free(BY);
        gsl_vector_free(BZ);

        //gsl_vector_free(bp);

        gsl_matrix_free(covX);
        gsl_matrix_free(covY);
        gsl_matrix_free(covZ);

        return lines;

    }

    static std::vector<Z3i::RealPoint> splines(std::vector<Z3i::RealPoint> points, double res){
        int nbPointForSample = points.size();
        double* T = (double*)malloc(nbPointForSample*sizeof(double));
        double* X = (double*)malloc(nbPointForSample*sizeof(double));
        double* Y = (double*)malloc(nbPointForSample*sizeof(double));
        double* Z = (double*)malloc(nbPointForSample*sizeof(double));
        //splines after remove point with detours
        for(int i = 0; i< nbPointForSample; i++){
            T[i] = i;
            X[i] = points[i][0];
            Y[i] = points[i][1];
            Z[i] = points[i][2];
        }

        gsl_interp_accel *accX = gsl_interp_accel_alloc ();
        gsl_spline *splineX = gsl_spline_alloc (gsl_interp_cspline, nbPointForSample);

        gsl_interp_accel *accY = gsl_interp_accel_alloc ();
        gsl_spline *splineY = gsl_spline_alloc (gsl_interp_cspline, nbPointForSample);

        gsl_interp_accel *accZ = gsl_interp_accel_alloc ();
        gsl_spline *splineZ = gsl_spline_alloc (gsl_interp_cspline, nbPointForSample);

        gsl_spline_init (splineX, T, X, nbPointForSample);
        gsl_spline_init (splineY, T, Y, nbPointForSample);
        gsl_spline_init (splineZ, T, Z, nbPointForSample);

        std::vector<Z3i::RealPoint> smoothFib;
        for(double t = 0; t <= nbPointForSample - 1; t += res){
            double x = gsl_spline_eval(splineX, t, accX);
            double y = gsl_spline_eval(splineY, t, accY);
            double z = gsl_spline_eval(splineZ, t, accZ);
            Z3i::RealPoint sPoint(x, y, z);
            smoothFib.push_back(sPoint);
        }

        gsl_spline_free (splineX);
        gsl_interp_accel_free (accX);

        gsl_spline_free (splineY);
        gsl_interp_accel_free (accY);

        gsl_spline_free (splineZ);
        gsl_interp_accel_free (accZ);

        free(T);
        free(X);
        free(Y);
        free(Z);        
        return smoothFib;
    }
private:
    static void bsplines(gsl_vector *xs, gsl_vector *ts, gsl_vector *wx, gsl_bspline_workspace *bw, gsl_vector *B, gsl_vector *c, gsl_matrix *cov, int nbCoeffs, int n) {
        //size_t nbBreak = nbCoeffs - 2;
        gsl_matrix *T;
        gsl_multifit_linear_workspace *mw;
        double chisq;

        gsl_vector *w;
        //w = gsl_vector_alloc(n);
        T = gsl_matrix_alloc(n, nbCoeffs);
        mw = gsl_multifit_linear_alloc(n, nbCoeffs);
        //construct the fit matrix X
        for (int i = 0; i < n; ++i){
            double ti = gsl_vector_get(ts, i);

            /* compute B_j(xi) for all j */
            gsl_bspline_eval(ti, B, bw);

            /* fill in row i of X */
            for (int j = 0; j < nbCoeffs; ++j){
                double Bj = gsl_vector_get(B, j);
                gsl_matrix_set(T, i, j, Bj);
            }
        }
        //gsl_multifit_linear(T, xs, c, cov, &chisq, mw);
        gsl_multifit_wlinear(T, wx, xs, c, cov, &chisq, mw);

        double Rsq, dof, tss;
        dof = n - nbCoeffs;
        tss = gsl_stats_wtss(wx->data, 1, xs->data, 1, ts->size);
        Rsq = 1.0 - chisq / tss;

        fprintf(stderr, "chisq/dof = %e, Rsq = %f\n", 
                chisq / dof, Rsq);

        //free 
        gsl_multifit_linear_free(mw);
        gsl_matrix_free(T);
        //gsl_vector_free(w);
    }
};
#endif
