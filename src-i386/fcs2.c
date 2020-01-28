// fcs2.c : C code for fast calculation of joint EQRs
//

#include <R.h>
#include <Rmath.h>
#include <R_ext/Utils.h>


// cdf
double pzinbinom(int q, double size, double nbmean, double zeroprob) {
	return zeroprob + (1.0 - zeroprob) * pnbinom_mu(q, size, nbmean, 1, 0);
}

// quantile function
int qzinbinom(double p, double size, double nbmean, double zeroprob) {
    // if zeroprob is 1, return 0
    if (zeroprob == 1.0)
        return 0;

    // calculate NB probability
    double prob;
    if (zeroprob >= p)
        prob = 0;
    else
        prob = (p - zeroprob) / (1 - zeroprob);

    // return NB quantile
    return qnbinom_mu(prob, size, nbmean, 1, 0);
}

// random sample
int rzinbinom(double size, double nbmean, double zeroprob) {
	if (unif_rand() < zeroprob)
		return 0;
	else
		return rnbinom_mu(size, nbmean);
}

// random constrained sample
int rzinbinom_constrained(double size, double nbmean, double zeroprob, int min, int max) {
    return qzinbinom(runif(pzinbinom(min - 1, size, nbmean, zeroprob),
                           pzinbinom(max, size, nbmean, zeroprob)),
                     size, nbmean, zeroprob);
}



// Joint EQR
void jointEQR(int catchMin[], int catchMax[], double size[], double nbmean[], double zeroprob[],
              int *nSurveys, int *nFits, int *nSamples, int *nSims, int *showProgress, double ret[])
{
	// Fetch random seed
	GetRNGstate();

	// declare variables
	int sample;			    // a MC sample
	double obsJointProb;	// observed joint probability
	double sampleJointProb; // sample joint probability
	int nSamplesBelowObs;	// counter for number of sample joint probs below observed joint prob

    int percentComplete = 0;  // percentage complete, for printing progress
    int lastPercentComplete = percentComplete;  // percentage complete, for printing progress

    // print progress
    if (*showProgress == 1) {
        Rprintf("0%%");
        R_FlushConsole();
    }

	for (int i = 0; i < *nSamples; ++i) {  // joint EQR sample i
		for (int j = 0; j < *nSurveys; ++j) {  // survey j
			// calculate observed joint probability:

			// calculate probability, for each fit
			// multiply over fits to calculate obs joint prob
			obsJointProb = 1.0;	 // reset to 1
			for (int l = 0; l < *nFits; ++l) {  // fit l
                // if catchMin < catchMax, sample catch
                if (catchMin[j + l * (*nSurveys)] < catchMax[j + l * (*nSurveys)]) {
                    sample = rzinbinom_constrained(size[i + l * (*nSamples)],
                                                   nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                   zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                   catchMin[j + l * (*nSurveys)],
                                                   catchMax[j + l * (*nSurveys)]);
                    obsJointProb *= pzinbinom(sample,
                                              size[i + l * (*nSamples)],
                                              nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                              zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)]);

                } else {
                    obsJointProb *= pzinbinom(catchMin[j + l * (*nSurveys)],
                                              size[i + l * (*nSamples)],
                                              nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                              zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                }
            }

            // calculate proportion of sampled join probabilities below observed:

            nSamplesBelowObs = 0;  // reset # samples below observed joint probability
			for (int k = 0; k < *nSims; ++k) {	// inner MC sample k
				// sample a ZINB obs, for each fit
				// calculate probability, for each fit
				// multiply over fits to calculate joint prob
				sampleJointProb = 1.0;	 // reset to 1
				for (int l = 0; l < *nFits; ++l) {  // fit l
					sample = rzinbinom(size[i + l * (*nSamples)],
									   nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
									   zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
					sampleJointProb *= pzinbinom(sample,
											     size[i + l * (*nSamples)],
									   			 nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
									   			 zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
				}

				// if sample <= observed joint prob, increase count
                if (sampleJointProb <= obsJointProb)  //  + 1e-10)
					++nSamplesBelowObs;
			}

			// result is then proportion of samples <= observed joint prob
			//   which equals # samples <= observed joint prob / # samples
            ret[i + j * (*nSamples)] = (double)nSamplesBelowObs / (double)(*nSims);  // ret[i, j]

            // check for user interruption
            R_CheckUserInterrupt();

            // update progress
            if (*showProgress == 1) {
                percentComplete = (100 * (j + i * (*nSurveys) + 1)) / ((*nSamples) * (*nSurveys));
                if (percentComplete > lastPercentComplete) {
                    lastPercentComplete = percentComplete;
                    Rprintf("\r%d%%", lastPercentComplete);
                    R_FlushConsole();
                }
            }
		}  // for survey j
	}  // for joint EQR sample i

    // end line of progress
    if (*showProgress == 1) {
        Rprintf("\n");
        R_FlushConsole();
    }

	// Put back random seed
	PutRNGstate();
}


// Joint EQR joining surveys as well as fits
void jointEQRJoiningSurveys(int catchMin[], int catchMax[], double size[], double nbmean[], double zeroprob[],
                            int iJoin[], int *nSurveys, int *nFits, int *nSamples, int *nSims, int *nJoinLevels, int *showProgress, double ret[])
{
    // Fetch random seed
    GetRNGstate();

    // declare variables
    int sample;             // a MC sample
    double obsJointProb;    // observed joint probability
    double sampleJointProb; // sample joint probability
    int nSamplesBelowObs;   // counter for number of sample joint probs below observed joint prob
    int nSurveysInSet;      // counter for number of surveys with current join level m

    // allocate memory to store surveys with join level m (of size *nSurveys as this is max)
    int *surveysInSet = (int *)R_alloc(*nSurveys, sizeof(int));

    int percentComplete = 0;  // percentage complete, for printing progress
    int lastPercentComplete = percentComplete;  // percentage complete, for printing progress

    // print progress
    if (*showProgress == 1) {
        Rprintf("0%%");
        R_FlushConsole();
    }

    for (int m = 0; m < *nJoinLevels; ++m) {  // join level (eg water body) m

        // create vector of surveys with join level m
        nSurveysInSet = 0;  // reset to 0
        for (int j = 0; j < *nSurveys; ++j)  // survey j
            if (iJoin[j] - 1 == m) {
                surveysInSet[nSurveysInSet] = j;
                ++nSurveysInSet;
            }

        // if no surveys, set all MC samples to NA
        if (nSurveysInSet == 0) {
            for (int i = 0; i < *nSamples; ++i)  // joint EQR sample i
                ret[i + m * (*nSamples)] = NA_REAL;

            // check for user interruption
            R_CheckUserInterrupt();

            // update progress
            if (*showProgress == 1) {
                percentComplete = (100 * (m + 1)) / (*nJoinLevels);
                if (percentComplete > lastPercentComplete) {
                    lastPercentComplete = percentComplete;
                    Rprintf("\r%d%%", lastPercentComplete);
                    R_FlushConsole();
                }
            }

        } else {

            // else, loop over MC samples
            for (int i = 0; i < *nSamples; ++i) {  // joint EQR sample i

                // calculate observed joint probability:

                // calculate probability, for each fit and each survey in set
                // multiply over fits and surveys to calculate obs joint prob
                obsJointProb = 1.0;  // reset to 1
                for (int j = 0; j < nSurveysInSet; ++j)  // survey indicator j
                    for (int l = 0; l < *nFits; ++l) {  // fit l
                        // if catchMin < catchMax, sample catch
                        if (catchMin[surveysInSet[j] + l * (*nSurveys)] < catchMax[surveysInSet[j] + l * (*nSurveys)]) {
                            sample = rzinbinom_constrained(size[i + l * (*nSamples)],
                                                           nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                           zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                           catchMin[surveysInSet[j] + l * (*nSurveys)],
                                                           catchMax[surveysInSet[j] + l * (*nSurveys)]);
                            obsJointProb *= pzinbinom(sample,
                                                      size[i + l * (*nSamples)],
                                                      nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                      zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);

                        } else {
                            obsJointProb *= pzinbinom(catchMin[surveysInSet[j] + l * (*nSurveys)],
                                                      size[i + l * (*nSamples)],
                                                      nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                      zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                        }
                    }

                // calculate proportion of sampled join probabilities below observed:

                nSamplesBelowObs = 0;  // reset # samples below observed joint probability
                for (int k = 0; k < *nSims; ++k) {  // inner MC sample k
                    // sample a ZINB obs, for each fit and each survey in set
                    // calculate probability, for each fit and each survey in set
                    // multiply over fits and surveys to calculate joint prob
                    sampleJointProb = 1.0;   // reset to 1
                    for (int j = 0; j < nSurveysInSet; ++j)  // survey indicator j
                        for (int l = 0; l < *nFits; ++l) {  // fit l
                            sample = rzinbinom(size[i + l * (*nSamples)],
                                               nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                               zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                            sampleJointProb *= pzinbinom(sample,
                                                         size[i + l * (*nSamples)],
                                                         nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                         zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                        }

                    // if sample <= observed joint prob, increase count
                    if (sampleJointProb <= obsJointProb)  //  + 1e-10)
                        ++nSamplesBelowObs;
                }

                // result is then proportion of samples <= observed joint prob
                //   which equals # samples <= observed joint prob / # samples
                ret[i + m * (*nSamples)] = (double)nSamplesBelowObs / (double)(*nSims);  // ret[i, m]

                // check for user interruption
                R_CheckUserInterrupt();

                // update progress
                if (*showProgress == 1) {
                    percentComplete = (100 * (i + m * (*nSamples) + 1)) / ((*nSamples) * (*nJoinLevels));
                    if (percentComplete > lastPercentComplete) {
                        lastPercentComplete = percentComplete;
                        Rprintf("\r%d%%", lastPercentComplete);
                        R_FlushConsole();
                    }
                }

            }  // for each joint EQR sample i
        }
    }  // for each join level m

    // end line of progress
    if (*showProgress == 1) {
        Rprintf("\n");
        R_FlushConsole();
    }

    // Put back random seed
    PutRNGstate();
}


// Joint EQR joining surveys and not
void jointEQRJoiningSurveysAndNot(int catchMin[], int catchMax[], double size[], double nbmean[], double zeroprob[],
                                  int iJoin[], int *nSurveys, int *nFits, int *nSamples, int *nSims, int *nJoinLevels,
                                  int *showProgress, double retBySurvey[], double retJoiningSurveys[])
{
    // Fetch random seed
    GetRNGstate();

    // declare variables
    int sample;             // a MC sample
    double obsJointProb;    // observed joint probability (over fits and surveys in set)
    double sampleJointProb; // sample joint probability
    int nSamplesBelowObs;   // counter for number of sample joint probs below observed joint prob
    int nSurveysInSet;      // counter for number of surveys with current join level m

    // allocate memory to store surveys with join level m (of size *nSurveys as this is max)
    int *surveysInSet = (int *)R_alloc(*nSurveys, sizeof(int));

    // allocate memory to store joint probabilities and sample count for surveys in current set
    double *obsJointProbPerSurvey = (double *)R_alloc(*nSurveys, sizeof(double));
    double *sampleJointProbPerSurvey = (double *)R_alloc(*nSurveys, sizeof(double));
    int *nSamplesBelowObsPerSurvey = (int *)R_alloc(*nSurveys, sizeof(int));

    int percentComplete = 0;  // percentage complete, for printing progress
    int lastPercentComplete = percentComplete;  // percentage complete, for printing progress

    // print progress
    if (*showProgress == 1) {
        Rprintf("0%%");
        R_FlushConsole();
    }

    // loop over join levels
    for (int m = 0; m < *nJoinLevels; ++m) {  // join level (eg water body) m

        // create vector of surveys with join level m
        nSurveysInSet = 0;  // reset to 0
        for (int j = 0; j < *nSurveys; ++j)  // survey j
            if (iJoin[j] - 1 == m) {
                surveysInSet[nSurveysInSet] = j;
                ++nSurveysInSet;
            }

        // if no surveys, set all MC samples to NA
        if (nSurveysInSet == 0) {
            for (int i = 0; i < *nSamples; ++i)  // joint EQR sample i
                retJoiningSurveys[i + m * (*nSamples)] = NA_REAL;

            // check for user interruption
            R_CheckUserInterrupt();

            // update progress
            if (*showProgress == 1) {
                percentComplete = (100 * (m + 1)) / (*nJoinLevels);
                if (percentComplete > lastPercentComplete) {
                    lastPercentComplete = percentComplete;
                    Rprintf("\r%d%%", lastPercentComplete);
                    R_FlushConsole();
                }
            }

        } else {

            // else, loop over MC samples
            for (int i = 0; i < *nSamples; ++i) {  // joint EQR sample i

                // calculate observed joint probability:

                // calculate probability, for each fit and each survey in set
                // multiply over fits and surveys to calculate obs joint prob
                obsJointProb = 1.0;  // reset to 1
                for (int j = 0; j < nSurveysInSet; ++j) {   // survey indicator j
                    obsJointProbPerSurvey[j] = 1.0;  // reset to 1

                    for (int l = 0; l < *nFits; ++l) {  // fit l
                        // if catchMin < catchMax, sample catch
                        if (catchMin[surveysInSet[j] + l * (*nSurveys)] < catchMax[surveysInSet[j] + l * (*nSurveys)]) {
                            sample = rzinbinom_constrained(size[i + l * (*nSamples)],
                                                           nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                           zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                           catchMin[surveysInSet[j] + l * (*nSurveys)],
                                                           catchMax[surveysInSet[j] + l * (*nSurveys)]);
                            obsJointProbPerSurvey[j] *= pzinbinom(sample,
                                                                  size[i + l * (*nSamples)],
                                                                  nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                                  zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);

                        } else {
                            obsJointProbPerSurvey[j] *= pzinbinom(catchMin[surveysInSet[j] + l * (*nSurveys)],
                                                                  size[i + l * (*nSamples)],
                                                                  nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                                  zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                        }
                    }

                    obsJointProb *= obsJointProbPerSurvey[j];
                }


                // calculate proportion of sampled join probabilities below observed:

                nSamplesBelowObs = 0;  // reset # samples below observed joint probability to 0
                for (int j = 0; j < nSurveysInSet; ++j)
                    nSamplesBelowObsPerSurvey[j]= 0;  // reset to 0
                for (int k = 0; k < *nSims; ++k) {  // inner MC sample k
                    // sample a ZINB obs, for each fit and each survey in set
                    // calculate probability, for each fit and each survey in set
                    // multiply over fits and surveys to calculate joint prob
                    sampleJointProb = 1.0;   // reset to 1
                    for (int j = 0; j < nSurveysInSet; ++j) {  // survey indicator j
                        sampleJointProbPerSurvey[j] = 1.0;  // reset to 1

                        for (int l = 0; l < *nFits; ++l) {  // fit l
                            sample = rzinbinom(size[i + l * (*nSamples)],
                                               nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                               zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                            sampleJointProbPerSurvey[j] *= pzinbinom(sample,
                                                                     size[i + l * (*nSamples)],
                                                                     nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                                     zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                        }

                        // if sample <= observed joint prob, increase count for this survey
                        if (sampleJointProbPerSurvey[j] <= obsJointProbPerSurvey[j])  // + 1.0e-10)
                            ++nSamplesBelowObsPerSurvey[j];

                        sampleJointProb *= sampleJointProbPerSurvey[j];
                    }

                    // if sample <= observed joint prob, increase count
                    if (sampleJointProb <= obsJointProb)  // + 1.0e-10)
                        ++nSamplesBelowObs;
                }

                // result is then proportion of samples <= observed joint prob
                //   which equals # samples <= observed joint prob / # samples
                retJoiningSurveys[i + m * (*nSamples)] = (double)nSamplesBelowObs / (double)(*nSims);  // retJoiningSurveys[i, m]
                for (int j = 0; j < nSurveysInSet; ++j)
                    retBySurvey[i + surveysInSet[j] * (*nSamples)] = (double)nSamplesBelowObsPerSurvey[j] / (double)(*nSims);

                // check for user interruption
                R_CheckUserInterrupt();

                // update progress
                if (*showProgress == 1) {
                    percentComplete = (100 * (i + m * (*nSamples) + 1)) / ((*nSamples) * (*nJoinLevels));
                    if (percentComplete > lastPercentComplete) {
                        lastPercentComplete = percentComplete;
                        Rprintf("\r%d%%", lastPercentComplete);
                        R_FlushConsole();
                    }
                }

            }  // for each joint EQR sample i
        }
    }  // for each join level m

    // end line of progress
    if (*showProgress == 1) {
        Rprintf("\n");
        R_FlushConsole();
    }

    // Put back random seed
    PutRNGstate();
}


// Joint and single EQRs
// NOTE: unlike joint EQR above, this function must check for NAs in catch or in nbmean & zeroprob
void jointAndSingleEQR(int catchMin[], int catchMax[], double size[], double nbmean[], double zeroprob[],
                       int *nSurveys, int *nFits, int *nSamples, int *nSims, int *showProgress, double ret[])
{
    // Fetch random seed
    GetRNGstate();

    // declare variables
    int sample;             // a MC sample
    double obsJointProb;    // observed joint probability
    double sampleJointProb; // sample joint probability
    int nSamplesBelowObs;   // counter for number of sample joint probs below observed joint prob

    // allocate memory to store probabilities and sample counts before joined
    double *obsProb = (double *)R_alloc(*nFits, sizeof(double));
    double *sampleProb = (double *)R_alloc(*nFits, sizeof(double));
    int *nSamplesBelowObsPerFit = (int *)R_alloc(*nFits, sizeof(int));

    int percentComplete = 0;  // percentage complete, for printing progress
    int lastPercentComplete = percentComplete;  // percentage complete, for printing progress

    // print progress
    if (*showProgress == 1) {
        Rprintf("0%%");
        R_FlushConsole();
    }

    for (int i = 0; i < *nSamples; ++i) {  // joint EQR sample i
        for (int j = 0; j < *nSurveys; ++j) {  // survey j
            // calculate observed joint probability:

            // calculate probability, for each fit
            // multiply over fits to calculate obs joint prob
            obsJointProb = 1.0;  // reset to 1
            for (int l = 0; l < *nFits; ++l) {  // fit l
                // check whether NA
                if (catchMin[j + l * (*nSurveys)] == NA_INTEGER ||
                    catchMax[j + l * (*nSurveys)] == NA_INTEGER ||
                    ISNA(nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)]) ||
                    ISNA(zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)])) {
                    // obs or prediction NA so set obs prob and joint prob to NA
                    obsProb[l] = obsJointProb = NA_REAL;

                } else {
                    // not NA so calculate
                    if (catchMin[j + l * (*nSurveys)] < catchMax[j + l * (*nSurveys)]) {
                        // catch min < max so sample
                        sample = rzinbinom_constrained(size[i + l * (*nSamples)],
                                                       nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                       zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                       catchMin[j + l * (*nSurveys)],
                                                       catchMax[j + l * (*nSurveys)]);
                        obsProb[l] = pzinbinom(sample,
                                               size[i + l * (*nSamples)],
                                               nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                               zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)]);

                    } else {
                        // catch min = max so calculate
                        obsProb[l] = pzinbinom(catchMin[j + l * (*nSurveys)],
                                               size[i + l * (*nSamples)],
                                               nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                               zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                    }

                    // update joint prob unless NA
                    if (!ISNA(obsJointProb))
                        obsJointProb *= obsProb[l];
                }
            }


            // calculate proportion of sampled join probabilities below observed:

            // check whether can calculate joint EQR
            if (ISNA(obsJointProb)) {
                // we cannot calculate joint EQR, so calculate singles where possible or set to NA

                // loop over fits
                for (int l = 0; l < *nFits; ++l) {  // fit l
                    // check whether can calculate single EQR
                    if (ISNA(obsProb[l])) {
                        // cannot calculate, so set EQR to NA
                        ret[i + j * (*nSamples) + (l + 1) * (*nSamples) * (*nSurveys)] = NA_REAL;

                    } else {
                        // can calculate

                        // reset # samples below observed prob
                        nSamplesBelowObsPerFit[l] = 0;  // reset to 0

                        for (int k = 0; k < *nSims; ++k) {  // inner MC sample k
                            sample = rzinbinom(size[i + l * (*nSamples)],
                                               nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                               zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                            sampleProb[l] = pzinbinom(sample,
                                                      size[i + l * (*nSamples)],
                                                      nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                      zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)]);

                            // if sample <= observed prob, increase count
                            if (sampleProb[l] <= obsProb[l])  //  + 1e-10)
                                ++nSamplesBelowObsPerFit[l];
                        }

                        // find single EQR from proportion
                        ret[i + j * (*nSamples) + (l + 1) * (*nSamples) * (*nSurveys)] = (double)nSamplesBelowObsPerFit[l] / (double)(*nSims);  // ret[i, j, l + 1] for fit l
                    }
                }

                // set joint EQR to NA
                ret[i + j * (*nSamples)] = NA_REAL;

            } else {
                // we can calculate joint EQR, so no need to check for NAs

                // reset # samples below observed joint probability
                nSamplesBelowObs = 0;  // reset to 0
                for (int l = 0; l < *nFits; ++l)
                    nSamplesBelowObsPerFit[l] = 0;  // reset to 0

                for (int k = 0; k < *nSims; ++k) {  // inner MC sample k
                    // sample a ZINB obs, for each fit
                    // calculate probability, for each fit
                    // multiply over fits to calculate joint prob
                    sampleJointProb = 1.0;   // reset to 1
                    for (int l = 0; l < *nFits; ++l) {  // fit l
                        sample = rzinbinom(size[i + l * (*nSamples)],
                                           nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                           zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                        sampleProb[l] = pzinbinom(sample,
                                                  size[i + l * (*nSamples)],
                                                  nbmean[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                  zeroprob[i + j * (*nSamples) + l * (*nSamples) * (*nSurveys)]);

                        // if sample <= observed prob, increase count
                        if (sampleProb[l] <= obsProb[l])  //  + 1e-10)
                            ++nSamplesBelowObsPerFit[l];

                        sampleJointProb *= sampleProb[l];
                    }

                    // if sample <= observed joint prob, increase count
                    if (sampleJointProb <= obsJointProb)  //  + 1e-10)
                        ++nSamplesBelowObs;
                }

                // result is then proportion of samples <= observed joint prob
                //   which equals # samples <= observed joint prob / # samples
                ret[i + j * (*nSamples)] = (double)nSamplesBelowObs / (double)(*nSims);  // ret[i, j, 1] for all fit EQR
                for (int l = 0; l < *nFits; ++l)
                    ret[i + j * (*nSamples) + (l + 1) * (*nSamples) * (*nSurveys)] = (double)nSamplesBelowObsPerFit[l] / (double)(*nSims);  // ret[i, j, l + 1] for fit l
            }

            // check for user interruption
            R_CheckUserInterrupt();

            // update progress
            if (*showProgress == 1) {
                percentComplete = (100 * (j + i * (*nSurveys) + 1)) / ((*nSamples) * (*nSurveys));
                if (percentComplete > lastPercentComplete) {
                    lastPercentComplete = percentComplete;
                    Rprintf("\r%d%%", lastPercentComplete);
                    R_FlushConsole();
                }
            }

        }  // for survey j
    }  // for joint EQR sample i

    // end line of progress
    if (*showProgress == 1) {
        Rprintf("\n");
        R_FlushConsole();
    }

    // Put back random seed
    PutRNGstate();
}


// Joint and single EQRs joining surveys as well as fits
void jointAndSingleEQRJoiningSurveys(int catchMin[], int catchMax[], double size[], double nbmean[], double zeroprob[], int iJoin[],
                                     int *nSurveys, int *nFits, int *nSamples, int *nSims, int *nJoinLevels, int *showProgress, double ret[])
{
    // Fetch random seed
    GetRNGstate();

    // declare variables
    int sample;             // a MC sample
    double obsJointProb;    // observed joint probability
    double sampleJointProb; // sample joint probability
    int nSamplesBelowObs;   // counter for number of sample joint probs below observed joint prob
    int nSurveysInSet;      // counter for number of surveys with current join level m

    // allocate memory to store surveys with join level m (of size *nSurveys as this is max)
    int *surveysInSet = (int *)R_alloc(*nSurveys, sizeof(int));

    // allocate memory to store probabilities and sample counts before joined by fit (but after joined by survey)
    double *obsJointProbPerFit = (double *)R_alloc(*nFits, sizeof(double));
    double *sampleJointProbPerFit = (double *)R_alloc(*nFits, sizeof(double));
    int *nSamplesBelowObsPerFit = (int *)R_alloc(*nFits, sizeof(int));

    double jointProbPerSurvey;   // single variable calculating joint prob per survey, for calculating overall joint prob
    double probPerSurveyAndFit;  // single variable storing observed prob

    int percentComplete = 0;  // percentage complete, for printing progress
    int lastPercentComplete = percentComplete;  // percentage complete, for printing progress

    // print progress
    if (*showProgress == 1) {
        Rprintf("0%%");
        R_FlushConsole();
    }

    for (int m = 0; m < *nJoinLevels; ++m) {  // join level (eg water body) m

        // create vector of surveys with join level m
        nSurveysInSet = 0;  // reset to 0
        for (int j = 0; j < *nSurveys; ++j)  // survey j
            if (iJoin[j] - 1 == m) {
                surveysInSet[nSurveysInSet] = j;
                ++nSurveysInSet;
            }

        // if no surveys, set all MC samples to NA
        if (nSurveysInSet == 0) {
            for (int i = 0; i < *nSamples; ++i)  // joint EQR sample i
                for (int l = 0; l <= *nFits; ++l)  // fit l (plus all fits)
                    ret[i + m * (*nSamples) + l * (*nSamples) * (*nJoinLevels)] = NA_REAL;

            // check for user interruption
            R_CheckUserInterrupt();

             // update progress
            if (*showProgress == 1) {
                percentComplete = (100 * (m + 1)) / (*nJoinLevels);
                if (percentComplete > lastPercentComplete) {
                    lastPercentComplete = percentComplete;
                    Rprintf("\r%d%%", lastPercentComplete);
                    R_FlushConsole();
                }
            }

        } else {

            // else, loop over MC samples
            for (int i = 0; i < *nSamples; ++i) {  // joint EQR sample i

                // calculate observed joint probability:

                // calculate probability, for each fit and each survey in set
                // multiply over fits and surveys to calculate obs joint prob
                obsJointProb = NA_REAL;  // set to NA initially
                for (int l = 0; l < *nFits; ++l)  // fit l
                    obsJointProbPerFit[l] = NA_REAL;  // set to NA initially
                for (int j = 0; j < nSurveysInSet; ++j) {  // survey indicator j
                    jointProbPerSurvey = 1.0;    // reset to 1

                    for (int l = 0; l < *nFits; ++l) {  // fit l
                        // check whether NA
                        if (catchMin[surveysInSet[j] + l * (*nSurveys)] == NA_INTEGER ||
                            catchMax[surveysInSet[j] + l * (*nSurveys)] == NA_INTEGER ||
                            ISNA(nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]) ||
                            ISNA(zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)])) {
                            // current or previous obs or prediction NA, so set joint prob per survey to NA
                            jointProbPerSurvey = NA_REAL;

                        } else {
                            // not NA so update
                            if (catchMin[surveysInSet[j] + l * (*nSurveys)] < catchMax[surveysInSet[j] + l * (*nSurveys)]) {
                                // catch min < max so sample
                                sample = rzinbinom_constrained(size[i + l * (*nSamples)],
                                                               nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                               zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                               catchMin[surveysInSet[j] + l * (*nSurveys)],
                                                               catchMax[surveysInSet[j] + l * (*nSurveys)]);
                                probPerSurveyAndFit = pzinbinom(sample,
                                                                size[i + l * (*nSamples)],
                                                                nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                                zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);

                            } else {
                                // catch min = max so calculate
                                probPerSurveyAndFit = pzinbinom(catchMin[surveysInSet[j] + l * (*nSurveys)],
                                                                size[i + l * (*nSamples)],
                                                                nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                                zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                            }

                            // update joint probs, unless NA
                            if (!ISNA(jointProbPerSurvey))
                                jointProbPerSurvey *= probPerSurveyAndFit;
                            if (ISNA(obsJointProbPerFit[l]))
                                obsJointProbPerFit[l] = probPerSurveyAndFit;
                            else
                                obsJointProbPerFit[l] *= probPerSurveyAndFit;
                        }
                    }

                    // update joint prob unless NA
                    if (!ISNA(jointProbPerSurvey)) {
                        if (ISNA(obsJointProb))
                            obsJointProb = jointProbPerSurvey;
                        else
                            obsJointProb *= jointProbPerSurvey;
                    }
                }


                // calculate proportion of sampled join probabilities below observed:

                // check whether can calculate joint EQR
                if (ISNA(obsJointProb)) {
                    // we cannot calculate joint EQR, so calculate singles where possible or set to NA

                    // loop over fits
                    for (int l = 0; l < *nFits; ++l) {  // fit l
                        // check whether can calculate single EQR
                        if (ISNA(obsJointProbPerFit[l])) {
                            // cannot calculate, so set EQR to NA
                            ret[i + m * (*nSamples) + (l + 1) * (*nSamples) * (*nJoinLevels)] = NA_REAL;

                        } else {
                            // can calculate, though not necessarily from all surveys

                            // reset # samples below observed joint probability
                            nSamplesBelowObsPerFit[l] = 0;  // reset to 0

                            for (int k = 0; k < *nSims; ++k) {  // inner MC sample k
                                sampleJointProbPerFit[l] = 1.0;  // reset to 1

                                for (int j = 0; j < nSurveysInSet; ++j) {  // survey indicator j
                                    // check whether NA (as not stored per survey and fit)
                                    if (catchMin[surveysInSet[j] + l * (*nSurveys)] != NA_INTEGER &&
                                        catchMax[surveysInSet[j] + l * (*nSurveys)] != NA_INTEGER &&
                                        !ISNA(nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]) &&
                                        !ISNA(zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)])) {
                                        sample = rzinbinom(size[i + l * (*nSamples)],
                                                           nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                           zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                                        sampleJointProbPerFit[l] *= pzinbinom(sample,
                                                                              size[i + l * (*nSamples)],
                                                                              nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                                              zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                                    }
                                }

                                // if sample <= observed prob, increase count
                                if (sampleJointProbPerFit[l] <= obsJointProbPerFit[l])  //  + 1e-10)
                                    ++nSamplesBelowObsPerFit[l];
                            }

                            // result is then proportion of samples <= observed joint prob
                            //   which equals # samples <= observed joint prob / # samples
                            ret[i + m * (*nSamples) + (l + 1) * (*nSamples) * (*nJoinLevels)] = (double)nSamplesBelowObsPerFit[l] / (double)(*nSims);  // ret[i, m, l + 1] for fit l
                        }
                    }  // fit l

                    // set all fits EQR to NA
                    ret[i + m * (*nSamples)] = NA_REAL;

                } else {
                    // we can calculate joint EQR, though not necessarily from all surveys

                    // reset # samples below observed joint probability
                    nSamplesBelowObs = 0;  // reset to 0
                    for (int l = 0; l < *nFits; ++l)
                        nSamplesBelowObsPerFit[l] = 0;  // reset to 0

                    for (int k = 0; k < *nSims; ++k) {  // inner MC sample k
                        // sample a ZINB obs, for each fit and each survey in set
                        // calculate probability, for each fit and each survey in set
                        // multiply over fits and surveys to calculate joint prob
                        sampleJointProb = 1.0;   // reset to 1
                        for (int l = 0; l < *nFits; ++l)  // fit l
                            sampleJointProbPerFit[l] = 1.0;  // reset to 1
                        for (int j = 0; j < nSurveysInSet; ++j) {  // survey indicator j
                            jointProbPerSurvey = 1.0;    // reset to 1

                            for (int l = 0; l < *nFits; ++l) {  // fit l
                                // check whether NA
                                if (catchMin[surveysInSet[j] + l * (*nSurveys)] == NA_INTEGER ||
                                    catchMax[surveysInSet[j] + l * (*nSurveys)] == NA_INTEGER ||
                                    ISNA(nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]) ||
                                    ISNA(zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)])) {
                                    // current or previous obs or prediction NA, so set joint prob per survey to NA
                                    jointProbPerSurvey = NA_REAL;

                                } else {
                                    sample = rzinbinom(size[i + l * (*nSamples)],
                                                       nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                       zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                                    probPerSurveyAndFit = pzinbinom(sample,
                                                                    size[i + l * (*nSamples)],
                                                                    nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                                    zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);

                                    // update joint probs, unless NA
                                    if (!ISNA(jointProbPerSurvey))
                                        jointProbPerSurvey *= probPerSurveyAndFit;
                                    if (!ISNA(obsJointProbPerFit[l]))
                                        sampleJointProbPerFit[l] *= probPerSurveyAndFit;
                                }
                            }

                            // update joint prob unless NA
                            if (!ISNA(jointProbPerSurvey))
                                sampleJointProb *= jointProbPerSurvey;
                        }

                        for (int l = 0; l < *nFits; ++l) {  // fit l
                            // if sample <= observed prob, increase count
                            if (!ISNA(obsJointProbPerFit[l]) && sampleJointProbPerFit[l] <= obsJointProbPerFit[l])  //  + 1e-10)
                                ++nSamplesBelowObsPerFit[l];
                        }

                        // if sample <= observed joint prob, increase count
                        if (sampleJointProb <= obsJointProb)  //  + 1e-10)
                            ++nSamplesBelowObs;
                    }

                    // result is then proportion of samples <= observed joint prob
                    //   which equals # samples <= observed joint prob / # samples
                    ret[i + m * (*nSamples)] = (double)nSamplesBelowObs / (double)(*nSims);  // ret[i, m, 1] for all fit EQR
                    for (int l = 0; l < *nFits; ++l) {  // fit l
                        if (ISNA(obsJointProbPerFit[l]))
                            ret[i + m * (*nSamples) + (l + 1) * (*nSamples) * (*nJoinLevels)] = NA_REAL;
                        else
                            ret[i + m * (*nSamples) + (l + 1) * (*nSamples) * (*nJoinLevels)] = (double)nSamplesBelowObsPerFit[l] / (double)(*nSims);  // ret[i, m, l + 1] for fit l
                    }
                }

                // check for user interruption
                R_CheckUserInterrupt();

                // update progress
                if (*showProgress == 1) {
                    percentComplete = (100 * (i + m * (*nSamples) + 1)) / ((*nSamples) * (*nJoinLevels));
                    if (percentComplete > lastPercentComplete) {
                        lastPercentComplete = percentComplete;
                        Rprintf("\r%d%%", lastPercentComplete);
                        R_FlushConsole();
                    }
                }

            }  // for each joint EQR sample i
        }
    }  // for each join level m

    // end line of progress
    if (*showProgress == 1) {
        Rprintf("\n");
        R_FlushConsole();
    }

    // Put back random seed
    PutRNGstate();
}


// Joint and single EQRs joining surveys and not
void jointAndSingleEQRJoiningSurveysAndNot(int catchMin[], int catchMax[], double size[], double nbmean[], double zeroprob[], int iJoin[],
                                           int *nSurveys, int *nFits, int *nSamples, int *nSims, int *nJoinLevels,
                                           int *showProgress, double retBySurvey[], double retJoiningSurveys[])
{
    // Fetch random seed
    GetRNGstate();

    // declare variables
    int sample;             // a MC sample
    double obsJointProb;    // observed joint probability (over fits and surveys in set)
    double sampleJointProb; // sample joint probability
    int nSamplesBelowObs;   // counter for number of sample joint probs below observed joint prob
    int nSurveysInSet;      // counter for number of surveys with current join level m

    // allocate memory to store surveys with join level m (of size *nSurveys as this is max)
    int *surveysInSet = (int *)R_alloc(*nSurveys, sizeof(int));

    // allocate memory to store probabilities and sample counts before joined (for each survey in set and each fit)
    double *obsProbPerSurveyAndFit = (double *)R_alloc((*nSurveys) * (*nFits), sizeof(double));
    double *sampleProbPerSurveyAndFit = (double *)R_alloc((*nSurveys) * (*nFits), sizeof(double));
    int *nSamplesBelowObsPerSurveyAndFit = (int *)R_alloc((*nSurveys) * (*nFits), sizeof(int));

    // allocate memory to store joint probabilities and sample count for surveys in current set
    double *obsJointProbPerSurvey = (double *)R_alloc(*nSurveys, sizeof(double));
    double *sampleJointProbPerSurvey = (double *)R_alloc(*nSurveys, sizeof(double));
    int *nSamplesBelowObsPerSurvey = (int *)R_alloc(*nSurveys, sizeof(int));

    // allocate memory to store probabilities and sample counts before joined by fit (but after joined by survey)
    double *obsJointProbPerFit = (double *)R_alloc(*nFits, sizeof(double));
    double *sampleJointProbPerFit = (double *)R_alloc(*nFits, sizeof(double));
    int *nSamplesBelowObsPerFit = (int *)R_alloc(*nFits, sizeof(int));

    int percentComplete = 0;  // percentage complete, for printing progress
    int lastPercentComplete = percentComplete;  // percentage complete, for printing progress

    // print progress
    if (*showProgress == 1) {
        Rprintf("0%%");
        R_FlushConsole();
    }

    // loop over join levels
    for (int m = 0; m < *nJoinLevels; ++m) {  // join level (eg water body) m

        // create vector of surveys with join level m
        nSurveysInSet = 0;  // reset to 0
        for (int j = 0; j < *nSurveys; ++j)  // survey j
            if (iJoin[j] - 1 == m) {
                surveysInSet[nSurveysInSet] = j;
                ++nSurveysInSet;
            }

        // if no surveys, set all MC samples to NA
        if (nSurveysInSet == 0) {
            for (int i = 0; i < *nSamples; ++i)  // joint EQR sample i
                for (int l = 0; l <= *nFits; ++l)  // fit l (plus all fits)
                    retJoiningSurveys[i + m * (*nSamples) + l * (*nSamples) * (*nJoinLevels)] = NA_REAL;

            // check for user interruption
            R_CheckUserInterrupt();

            // update progress
            if (*showProgress == 1) {
                percentComplete = (100 * (m + 1)) / (*nJoinLevels);
                if (percentComplete > lastPercentComplete) {
                    lastPercentComplete = percentComplete;
                    Rprintf("\r%d%%", lastPercentComplete);
                    R_FlushConsole();
                }
            }

        } else {

            // else, loop over MC samples
            for (int i = 0; i < *nSamples; ++i) {  // joint EQR sample i

                // calculate observed joint probability:

                // calculate probability, for each fit and each survey in set
                // multiply over fits and surveys to calculate obs joint prob

                // reset joint probs to 1
                obsJointProb = NA_REAL;  // set to NA initially
                for (int l = 0; l < *nFits; ++l)  // fit l
                    obsJointProbPerFit[l] = NA_REAL;  // set to NA initially
                for (int j = 0; j < nSurveysInSet; ++j) {  // survey indicator j
                    obsJointProbPerSurvey[j] = 1.0;  // reset to 1

                    for (int l = 0; l < *nFits; ++l) {  // fit l
                        // check whether NA
                        if (catchMin[surveysInSet[j] + l * (*nSurveys)] == NA_INTEGER ||
                            catchMax[surveysInSet[j] + l * (*nSurveys)] == NA_INTEGER ||
                            ISNA(nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]) ||
                            ISNA(zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)])) {
                            // current obs or prediction NA, so set obs and joint prob per survey to NA
                            obsProbPerSurveyAndFit[j + l * nSurveysInSet] = obsJointProbPerSurvey[j] = NA_REAL;

                        } else {
                            // not NA so calculate
                            if (catchMin[surveysInSet[j] + l * (*nSurveys)] < catchMax[surveysInSet[j] + l * (*nSurveys)]) {
                                // catch min < max so sample
                                sample = rzinbinom_constrained(size[i + l * (*nSamples)],
                                                               nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                               zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                               catchMin[surveysInSet[j] + l * (*nSurveys)],
                                                               catchMax[surveysInSet[j] + l * (*nSurveys)]);
                                obsProbPerSurveyAndFit[j + l * nSurveysInSet] = pzinbinom(sample,
                                                                                          size[i + l * (*nSamples)],
                                                                                          nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                                                          zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);

                            } else {
                                // catch min = max so calculate
                                obsProbPerSurveyAndFit[j + l * nSurveysInSet] = pzinbinom(catchMin[surveysInSet[j] + l * (*nSurveys)],
                                                                                          size[i + l * (*nSamples)],
                                                                                          nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                                                          zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                            }

                            // update joint probs, unless NA
                            if (!ISNA(obsJointProbPerSurvey[j]))
                                obsJointProbPerSurvey[j] *= obsProbPerSurveyAndFit[j + l * nSurveysInSet];
                            if (ISNA(obsJointProbPerFit[l]))
                                obsJointProbPerFit[l] = obsProbPerSurveyAndFit[j + l * nSurveysInSet];
                            else
                                obsJointProbPerFit[l] *= obsProbPerSurveyAndFit[j + l * nSurveysInSet];
                        }
                    }

                    // update overall joint prob, unless NA
                    if (!ISNA(obsJointProbPerSurvey[j])) {
                        if (ISNA(obsJointProb))
                            obsJointProb = obsJointProbPerSurvey[j];
                        else
                            obsJointProb *= obsJointProbPerSurvey[j];
                    }
                }


                // calculate proportion of sampled join probabilities below observed:

                // reset # samples below observed joint probability to 0
                nSamplesBelowObs = 0;
                for (int j = 0; j < nSurveysInSet; ++j) { // survey indicator j
                    nSamplesBelowObsPerSurvey[j]= 0;
                    for (int l = 0; l < *nFits; ++l) {  // fit l
                        nSamplesBelowObsPerSurveyAndFit[j + l * nSurveysInSet] = 0;
                        nSamplesBelowObsPerFit[l]= 0;
                    }
                }

                // loop over MC samples
                for (int k = 0; k < *nSims; ++k) {  // inner MC sample k
                    // sample a ZINB obs, for each fit and each survey in set
                    // calculate probability, for each fit and each survey in set
                    // multiply over fits and surveys to calculate joint prob

                    // reset joint probs to 1
                    sampleJointProb = 1.0;   // reset to 1
                    for (int l = 0; l < *nFits; ++l)  // fit l
                        sampleJointProbPerFit[l] = 1.0;   // reset to 1
                    for (int j = 0; j < nSurveysInSet; ++j) {  // survey indicator j
                        sampleJointProbPerSurvey[j] = 1.0;  // reset to 1

                        for (int l = 0; l < *nFits; ++l) {  // fit l
                            // check whether this survey & fit OK
                            if (!ISNA(obsProbPerSurveyAndFit[j + l * nSurveysInSet])) {
                                sample = rzinbinom(size[i + l * (*nSamples)],
                                                   nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                   zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);
                                sampleProbPerSurveyAndFit[j + l * nSurveysInSet] = pzinbinom(sample,
                                                                                             size[i + l * (*nSamples)],
                                                                                             nbmean[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)],
                                                                                             zeroprob[i + surveysInSet[j] * (*nSamples) + l * (*nSamples) * (*nSurveys)]);

                                // update joint probs, if possible
                                if (!ISNA(obsJointProbPerFit[l]))
                                    sampleJointProbPerFit[l] *= sampleProbPerSurveyAndFit[j + l * nSurveysInSet];
                                if (!ISNA(obsJointProbPerSurvey[j]))
                                    sampleJointProbPerSurvey[j] *= sampleProbPerSurveyAndFit[j + l * nSurveysInSet];

                                // if sample <= observed prob, increase count for this survey and fit
                                if (sampleProbPerSurveyAndFit[j + l * nSurveysInSet] <= obsProbPerSurveyAndFit[j + l * nSurveysInSet])  // + 1.0e-10)
                                    ++nSamplesBelowObsPerSurveyAndFit[j + l * nSurveysInSet];
                            }
                        }

                        // if sample <= observed joint prob, increase count for this survey
                        if (!ISNA(obsJointProbPerSurvey[j]) && sampleJointProbPerSurvey[j] <= obsJointProbPerSurvey[j])  // + 1.0e-10)
                            ++nSamplesBelowObsPerSurvey[j];

                        // update overall joint prob, if possible
                        if (!ISNA(obsJointProb))
                            sampleJointProb *= sampleJointProbPerSurvey[j];
                    }

                    for (int l = 0; l < *nFits; ++l) {  // fit l
                        // if sample <= observed joint prob, increase count for this fit
                        if (!ISNA(obsJointProbPerFit[l]) && sampleJointProbPerFit[l] <= obsJointProbPerFit[l])  // + 1.0e-10)
                            ++nSamplesBelowObsPerFit[l];
                    }

                    // if sample <= observed joint prob, increase count
                    if (!ISNA(obsJointProb) && sampleJointProb <= obsJointProb)  // + 1.0e-10)
                        ++nSamplesBelowObs;
                } // for inner MC sample k

                // result is then proportion of samples <= observed joint prob
                //   which equals # samples <= observed joint prob / # samples
                if (ISNA(obsJointProb))
                    retJoiningSurveys[i + m * (*nSamples)] = NA_REAL;
                else
                    retJoiningSurveys[i + m * (*nSamples)] = (double)nSamplesBelowObs / (double)(*nSims);  // retJoiningSurveys[i, m, 1] for all fit EQR

                for (int j = 0; j < nSurveysInSet; ++j) {
                    if (ISNA(obsJointProbPerSurvey[j]))
                        retBySurvey[i + surveysInSet[j] * (*nSamples)] = NA_REAL;
                    else
                        retBySurvey[i + surveysInSet[j] * (*nSamples)] = (double)nSamplesBelowObsPerSurvey[j] / (double)(*nSims);  // retBySurvey[i, surveysInSet[j], 1] for all fit EQR
                }

                for (int l = 0; l < *nFits; ++l) {  // fit l
                    if (ISNA(obsJointProbPerFit[l]))
                        retJoiningSurveys[i + m * (*nSamples) + (l + 1) * (*nSamples) * (*nJoinLevels)] = NA_REAL;
                    else
                        retJoiningSurveys[i + m * (*nSamples) + (l + 1) * (*nSamples) * (*nJoinLevels)] =
                                    (double)nSamplesBelowObsPerFit[l] / (double)(*nSims);  // retJoiningSurveys[i, m, l + 1] for fit l

                    for (int j = 0; j < nSurveysInSet; ++j) {
                        if (ISNA(obsProbPerSurveyAndFit[j + l * nSurveysInSet]))
                            retBySurvey[i + surveysInSet[j] * (*nSamples) + (l + 1) * (*nSamples) * (*nSurveys)] = NA_REAL;
                        else
                            retBySurvey[i + surveysInSet[j] * (*nSamples) + (l + 1) * (*nSamples) * (*nSurveys)] =
                                    (double)nSamplesBelowObsPerSurveyAndFit[j + l * nSurveysInSet] / (double)(*nSims);  // retBySurvey[i, surveysInSet[j], l + 1] for fit l
                    }
                }

                // check for user interruption
                R_CheckUserInterrupt();

                // update progress
                if (*showProgress == 1) {
                    percentComplete = (100 * (i + m * (*nSamples) + 1)) / ((*nSamples) * (*nJoinLevels));
                    if (percentComplete > lastPercentComplete) {
                        lastPercentComplete = percentComplete;
                        Rprintf("\r%d%%", lastPercentComplete);
                        R_FlushConsole();
                    }
                }

            }  // for each joint EQR sample i
        }
    }  // for each join level m

    // end line of progress
    if (*showProgress == 1) {
        Rprintf("\n");
        R_FlushConsole();
    }

    // Put back random seed
    PutRNGstate();
}

