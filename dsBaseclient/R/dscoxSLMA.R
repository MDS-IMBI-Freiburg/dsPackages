#'
#' @title Fit a cox proportional hazard Model (coxph) with pooling via Study Level Meta-Analysis (SLMA)
#' @description Fits a cox proportional hazard Model (coxph) on data from single or multiple sources
#' with pooled co-analysis across studies being based on SLMA (Study Level Meta Analysis).
#' @details \code{ds.coxSLMA} specifies the structure of a cox proportional hazard Model (coxph)
#' to be fitted separately on each study or data source. Calls the serverside functions
#' coxSLMADS1 (aggregate) and coxSLMADS2 (aggregate).ds.coxSLMA sends
#' a command to every data source to fit the model required but each separate source
#' simply fits that model to completion (ie undertakes all iterations until
#' the model converges) and the estimates (regression coefficients) and their standard
#' errors from each source are sent back to the client and are then pooled using SLMA
#' via any approach the user wishes to implement. The ds.coxSLMA functions includes
#' a function  <mixmeta> which  pools the models
#' across studies using the mixmeta function (from the mixmeta package) using three
#' optimisation methods: random effects under maximum likelihood (ML); random effects
#' under restricted maximum likelihood (REML); or fixed effects (FE). But once
#' the estimates and standard errors are on the clientside, the user
#' can alternatively choose to use the metafor package in any way he/she wishes,
#' to pool the coefficients across studies or, indeed, to use another
#' meta-analysis package, or their own code.
#' In \code{formula} Most shortcut notation for formulas allowed under R's standard \code{coxph()}
#' function is also allowed by \code{ds.coxSLMA}.
#'
#' coxph can be fitted very simply using a formula such as:
#'
#' \deqn{y~a+b+c+d}
#'
#' which simply means fit a coxph with \code{y} as the outcome variable (Survival object
#' calculated separately using ds.Surv and stored in the server) and
#' \code{a}, \code{b}, \code{c} and \code{d} as covariates.
#' By default all such models also include an intercept (regression constant) term.
#' The \code{dataName} argument avoids you having to specify the name of the
#' data frame in front of each covariate in the formula.
#' For example, if the data frame is called \code{DataFrame} you
#' avoid having to write: \eqn{DataFrame$y~DataFrame$a+DataFrame$b+DataFrame$c+DataFrame$d}
#'
#' The \code{checks} argument verifies that the variables in the model are all defined (exist)
#' on the server-site at every study
#' and that they have the correct characteristics required to fit the model.
#' It is suggested to make \code{checks} argument TRUE only if an unexplained
#'  problem in the model fit is encountered because the running process takes several minutes.
#'
#' In \code{maxit} Logistic regression and Poisson regression
#' models can require many iterations, particularly if the starting value of the
#' regression constant is far away from its actual value that the GLM
#' is trying to estimate. In consequence we often set \code{maxit=30}
#' but depending on the nature of the models you wish to fit, you may wish
#' to be alerted much more quickly than this if there is a delay in convergence,
#' or you may wish to allow more iterations.
#'
#' Server functions called: \code{coxSLMADS1}, and \code{coxSLMADS2}.
#' @param formula an object of class formula describing
#' the model to be fitted. For more information see
#' \strong{Details}.
#' @param weights a character string specifying the name of a variable containing
#' prior regression weights for the fitting process. \code{ds.coxSLMA} does not allow a weights vector to be
#' written directly into the GLM formula.
#' @param combine.with.metafor logical. If TRUE the
#' estimates and standard errors for each regression coefficient are pooled across
#' studies using random-effects meta-analysis under maximum likelihood (ML),
#' restricted maximum likelihood (REML) or fixed-effects meta-analysis (FE). Default TRUE.
#' @param dataName a character string specifying the name of an (optional) data frame
#' that contains all of the variables in the coxph formula.
#' @param checks logical. If TRUE \code{ds.coxSLMA} checks the structural integrity
#' of the model. Default FALSE. For more information see \strong{Details}.
#' @param maxit a numeric scalar denoting the maximum number of iterations that
#' are permitted before \code{ds.coxSLMA} declares that the model has failed to converge.
#' For more information see \strong{Details}.
#' @param datasources a list of \code{\link{DSConnection-class}} objects obtained after login.
#' If the \code{datasources} argument is not specified
#' the default set of connections will be used: see \code{\link{datashield.connections_default}}.
#' @return The serverside aggregate functions \code{coxSLMADS1} and \code{coxSLMADS2} return
#' output to the clientside.
#' This is precisely the same as the coxph object that is usually created by a call to coxph() in native R and it
#' contains all the same elements (see help for coxph in native R). Because it is a serverside
#' object, no disclosure blocks apply. However, such disclosure blocks do apply to the information
#' passed to the clientside. In consequence, rather than containing all the components of a
#' standard coxph object in native R, the components of the coxph object that are returned by
#' \code{ds.coxSLMA} include: a mixture of non-disclosive elements of the coxph object
#' reported separately by study included in a list object called \code{output.summary}; and
#' a series of other list objects that represent inferences aggregated across studies.
#' @return the study specific items include:
#' @return \code{coefficients}: a matrix with 5 columns:
#'    \itemize{
#'    \item{First}{: the names of all of the regression parameters (coefficients) in the model}
#'    \item{second}{: the estimated values}
#'    \item{third}{: corresponding standard errors of the estimated values}
#'    \item{fourth}{: the ratio of estimate/standard error}
#'    \item{fifth}{: the p-value treating that as a standardised normal deviate}
#' }
#' @return \code{formula}: model formula, see description of formula as an input parameter (above).
#' @return \code{df.resid}: the residual degrees of freedom around the model.
#' @return \code{deviance.resid}: the residual deviance around the model.
#' @return \code{df.null}: the degrees of freedom around the null model (with just an intercept).
#' @return \code{dev.null}: the deviance around the null model (with just an intercept).
#' @return \code{CorrMatrix}: the correlation matrix of parameter estimates.
#' @return \code{VarCovMatrix}: the variance-covariance matrix of parameter estimates.
#' @return \code{weights}: the name of the vector (if any) holding regression weights.
#' @return \code{cov.scaled}: equivalent to \code{VarCovMatrix}.
#' @return \code{cov.unscaled}: equivalent to VarCovMatrix but assuming dispersion (scale)
#' parameter is 1.
#' @return \code{Nmissing}: the number of missing observations in the given study.
#' @return \code{Nvalid}: the number of valid (non-missing) observations in the given study.
#' @return \code{Ntotal}: the total number of observations in the given study
#'                        (\code{Nvalid} + \code{Nmissing}).
#' @return \code{data}: equivalent to input parameter \code{dataName} (above).
#' @return \code{dispersion}: the estimated dispersion parameter: deviance.resid/df.resid for
#' a gaussian family multiple regression model, 1.00 for logistic and poisson regression.
#' @return \code{call}:  summary of key elements of the call to fit the model.
#' @return \code{na.action}:  chosen method of dealing with missing values. This is
#' usually, \code{na.action = na.omit} - see help in native R.
#' @return \code{iter}: the number of iterations required to achieve convergence
#' of the glm model in each separate study.
#' @return Once the study-specific output has been returned, \code{ds.coxSLMA}
#' returns a series of lists relating to the aggregated inferences across studies.
#' These include the following:
#' @return \code{num.valid.studies}: the number of studies with valid output
#' included in the combined analysis
#' @author Sofack, Ghislain.(Based on ds.glmSLMA by Paul Burton for DataSHIELD Development Team)

#' @examples
#' \dontrun{
#'
#'  ## Version 6, for version 5 see Wiki
#'   # Connecting to the Opal servers
#'
#'   require('DSI')
#'   require('DSOpal')
#'   require('dsBaseClient')
#'
#'   # Example 1: Fitting GLM for survival analysis
#'   # For this analysis we need to load survival data from the server
#'
#'   builder <- DSI::newDSLoginBuilder()
#'   builder$append(server = "study1",
#'                  url = "http://192.168.56.100:8080/",
#'                  user = "administrator", password = "datashield_test&",
#'                  table = "SURVIVAL.EXPAND_NO_MISSING1", driver = "OpalDriver")
#'   builder$append(server = "study2",
#'                  url = "http://192.168.56.100:8080/",
#'                  user = "administrator", password = "datashield_test&",
#'                  table = "SURVIVAL.EXPAND_NO_MISSING2", driver = "OpalDriver")
#'   builder$append(server = "study3",
#'                  url = "http://192.168.56.100:8080/",
#'                  user = "administrator", password = "datashield_test&",
#'                  table = "SURVIVAL.EXPAND_NO_MISSING3", driver = "OpalDriver")
#'   logindata <- builder$build()
#'
#'   # Log onto the remote Opal training servers
#'   connections <- DSI::datashield.login(logins = logindata, assign = TRUE, symbol = "D")
#'
#'   # Create the serverside survival object
#'
#'   ds.Surv (x = "D$survtime",
#'            y = "D$cens",
#'            newobj = "Survobj"
#'            datasources = connections)
#'
#'
#'   # make sure that the outcome is numeric
#'   ds.asNumeric(x.name = "D$cens",
#'                newobj = "EVENT",
#'                datasources = connections)
#'
#'
#'
#'
#'   ds.coxSLMA(formula = Survobj ~ female * age.60,
#'          data = "D",
#'          weights = NULL,
#'          checks = FALSE,
#'          maxit = 20,
#'          datasources = connections)
#'
#'   # Clear the Datashield R sessions and logout
#'   datashield.logout(connections)
#'
#' }
#'
#' @export








ds.coxSLMA<-function(formula=NULL, weights=NULL,  combine.with.metafor=TRUE,dataName=NULL,
                     checks=FALSE, maxit=30, datasources=NULL) {

  # look for DS connections
  if(is.null(datasources)){
    datasources <- datashield.connections_find()
  }


  # verify that 'formula' was set
  if(is.null(formula)){
    stop(" Please provide a valid regression formula!", call.=FALSE)
  }


  # check if user gave offset or weights directly in formula, if so the argument 'offset' or 'weights'
  # to provide name of offset or weights variable
  if(sum(as.numeric(grepl('weights', formula, ignore.case=TRUE)))>0)
  {
    cat("\n\n WARNING: you may have specified an offset or regression weights")
    cat("\n as part of the model formula. In ds.glm (unlike the usual glm in R)")
    cat("\n you must specify an offset or weights separately from the formula")
    cat("\n using the offset or weights argument.\n\n")
  }

  formula <- stats::as.formula(formula)


  # if the argument 'dataName' is set, check that the data frame is defined (i.e. exists) on the server site
  if(!(is.null(dataName))){
    defined <- dsBaseClient:::isDefined(datasources, dataName)
  }


  #MOVE ITERATION COUNT BEFORE ASSIGNMENT OF beta.vect.next
  #Iterations need to be counted. Start off with the count at 0
  #and increment by 1 at each new iteration

  iteration.count<-0

  # number of 'valid' studies (those that passed the checks) and vector of beta values
  numstudies <- length(datasources)


  #ARBITRARY LENGTH FOR START BETAs AT THIS STAGE BUT IN LEGAL TRANSMISSION FORMAT ("0,0,0,0,0")
  beta.vect.next <- rep(0,5)
  beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")


  #IDENTIFY THE CORRECT DIMENSION FOR START BETAs VIA CALLING FIRST COMPONENT OF glmDS

  cally1 <- call('coxSLMADS1', formula, weights, dataName)

  study.summary.0 <- DSI::datashield.aggregate(datasources, cally1)


  at.least.one.study.data.error<-0
  at.least.one.study.valid<-0



  num.par.coxph<-NULL

  coef.names<-NULL

  for(hh in 1:numstudies) {
    if(study.summary.0[[hh]]$errorMessage!="No errors")
    {
      at.least.one.study.data.error<-1
    }else{
      at.least.one.study.valid<-1
      num.par.coxph<-study.summary.0[[hh]][[1]][[2]]
      coef.names<-study.summary.0[[hh]][[2]]
    }
  }


  y.invalid<-NULL
  Xpar.invalid<-NULL
  w.invalid<-NULL
  coxph.saturation.invalid<-NULL
  errorMessage<-NULL



  for(ss in 1:numstudies)
  {
    y.invalid<-c(y.invalid,study.summary.0[[ss]][[3]])
    Xpar.invalid<-rbind(Xpar.invalid,study.summary.0[[ss]][[4]])
    w.invalid<-c(w.invalid,study.summary.0[[ss]][[5]])
    coxph.saturation.invalid <-c(coxph.saturation.invalid,study.summary.0[[ss]][[6]])
    errorMessage<-c(errorMessage,study.summary.0[[ss]][[7]])
  }




  y.invalid<-as.matrix(y.invalid)
  sum.y.invalid<-sum(y.invalid)
  dimnames(y.invalid)<-list(names(datasources),"Y VECTOR")

  Xpar.invalid<-as.matrix(Xpar.invalid)
  sum.Xpar.invalid<-sum(Xpar.invalid)
  dimnames(Xpar.invalid)<-list(names(datasources),coef.names)

  w.invalid<-as.matrix(w.invalid)
  sum.w.invalid<-sum(w.invalid)
  dimnames(w.invalid)<-list(names(datasources),"WEIGHT VECTOR")


  coxph.saturation.invalid<-as.matrix(coxph.saturation.invalid)
  sum.coxph.saturation.invalid<-sum(coxph.saturation.invalid)
  dimnames(coxph.saturation.invalid)<-list(names(datasources),"MODEL OVERPARAMETERIZED")

  errorMessage<-as.matrix(errorMessage)
  dimnames(errorMessage)<-list(names(datasources),"ERROR MESSAGES")



  output.blocked.information.1<-"EVERY STUDY HAS DATA THAT COULD BE POTENTIALLY DISCLOSIVE UNDER THE CURRENT MODEL:"
  output.blocked.information.2<-"Any values of 1 in the following tables denote potential disclosure risks."
  output.blocked.information.3<-"Please use the argument <datasources> to include only valid studies."
  output.blocked.information.4<-"Errors by study are as follows:"


  #CASE 1 - NO STUDIES VALID
  if(!at.least.one.study.valid)
  {
    message("\n\nEVERY STUDY HAS DATA THAT COULD BE POTENTIALLY DISCLOSIVE UNDER THE CURRENT MODEL:\n",
            "Any values of 1 in the following tables denote potential disclosure risks.\n",
            "Errors by study are as follows:\n")
    #		print(as.matrix(y.invalid))
    #		print(as.matrix(Xpar.invalid))
    #		print(as.matrix(w.invalid))
    #		print(as.matrix(coxph.saturation.invalid))
    #		print(as.matrix(errorMessage))


    return(list(
      output.blocked.information.1,
      output.blocked.information.2,
      output.blocked.information.4,
      y.vector.error=y.invalid,
      X.matrix.error=Xpar.invalid,
      weight.vector.error=w.invalid,
      coxph.overparameterized=coxph.saturation.invalid,
      errorMessage=errorMessage
    ))
  }

  #CASE 2 - AT LEAST ONE STUDY VALID AND AT LEAST ONE INVALID
  if(at.least.one.study.data.error)
  {
    message("\n\nAT LEAST ONE STUDY HAS DATA THAT COULD BE POTENTIALLY DISCLOSIVE UNDER THE CURRENT MODEL:\n",
            "Any values of 1 in the following tables denote potential disclosure risks.\n",
            "No analytic results are returned for potentially disclosive studies and\n",
            "pooled co-estimates across studies are based only on the valid studies.\n",
            "You may also choose to exclude invalid studies from\n",
            "the whole analysis using the <datasources> argument.\n",
            "Errors by study are as follows:\n")
    #		print(as.matrix(y.invalid))
    #		print(as.matrix(Xpar.invalid))
    #		print(as.matrix(w.invalid))
    #		print(as.matrix(coxph.saturation.invalid))
    #		print(as.matrix(errorMessage))

  }



  beta.vect.next <- rep(0,num.par.coxph)
  beta.vect.temp <- paste0(as.character(beta.vect.next), collapse=",")


  #Provide arbitrary starting value for deviance to enable subsequent calculation of the
  #change in deviance between iterations
  dev.old<-9.99e+99

  #Convergence state needs to be monitored.
  converge.state<-FALSE

  #Define a convergence criterion. This value of epsilon corresponds to that used
  #by default for GLMs in R (see section S3 for details)
  epsilon<-1.0e-08

  f<-NULL


  #NOW CALL SECOND COMPONENT OF glmDS TO GENERATE SCORE VECTORS AND INFORMATION MATRICES

  cally2 <- call('coxSLMADS2', formula, weights, dataName)

  study.summary <- DSI::datashield.aggregate(datasources, cally2)




  numstudies<-length(datasources)

  study.include.in.analysis <- NULL
  study.with.errors<-NULL
  all.studies.valid<-1
  no.studies.valid<-1



  for(j in 1:numstudies)
  {
    ss1<-unlist(study.summary[[j]][[1]])
    if(is.numeric(ss1))
    {
      inv.diag.se<-1/sqrt(diag(study.summary[[j]]$cov.scaled))

      cor.matrix<-t(diag(inv.diag.se))%*%study.summary[[j]]$cov.scaled%*%(diag(inv.diag.se))
      study.summary[[j]]$VarCovMatrix<-study.summary[[j]]$cov.scaled
      study.summary[[j]]$CorrMatrix<-cor.matrix
      study.include.in.analysis<-c(study.include.in.analysis,j)
      no.studies.valid<-0
    }else{
      study.with.errors<-c(study.with.errors,j)
      all.studies.valid<-0
    }
  }



  #MAKE SURE THAT IF SOME STUDIES HAVE MORE PARAMETERS IN THE
  #FITTED coxph (eg BECAUSE OF ALIASING) THE FINAL RETURN MATRICES
  #HAVE ENOUGH ROWS TO FIT THE MAXIMUM LENGTH


  numcoefficients.max<-0



  #
  # for(g in 1:numstudies){
  #
  #   coefs = rbind(study.summary[[1]]$coefficients[,1], study.summary[[2]]$coefficients[,1], study.summary[[3]]$coefficients[,1])
  #
  #   covars = list( study.summary[[1]]$VarCovMatrix,study.summary[[2]]$VarCovMatrix, study.summary[[3]]$VarCovMatrix )
  #
  #   if(numstudies > 1){
  #
  #   mix <- mixmeta::mixmeta(coefs, S = covars, method = "fixed")
  #
  #   return(list(study.summary = study.summary, coefs= coefs, covars = covars, mix = summary(mix)))
  #
  #   }
  #   else {
  #     return(list(study.summary = study.summary, coefs = coefs, covars = covars))
  #   }
  #
  # }
  #



  for(g in range(1:numstudies)){

    coefs = study.summary[[g]]$coefficients[,1]

    covars = list( study.summary[[g]]$VarCovMatrix )

    return(list(datasources = datasources, study.summary = study.summary, coefs= coefs, covars = covars))

  }






  # for(g in study.include.in.analysis){
  #   if(length(study.summary[[g]]$coefficients[,1])>numcoefficients.max){
  #     numcoefficients.max<-length(study.summary[[g]]$coefficients[,1])
  #   }
  # }
  #

  numcoefficients<-numcoefficients.max

  betamatrix<-matrix(NA,nrow<-numcoefficients,ncol=numstudies)

  # sematrix<-matrix(NA,nrow<-numcoefficients,ncol=numstudies)


  for(k in study.include.in.analysis){
    betamatrix[,k]<-study.summary[[k]]$coefficients[,1]

    # sematrix[,k]<-study.summary[[k]]$coefficients[,2]
  }

  return( coefs)

  ################################################
  #ANNOTATE OUTPUT MATRICES WITH STUDY INDICATORS#
  ################################################

  # study.names.list<-NULL
  # betas.study.names.list<-NULL
  # ses.study.names.list<-NULL
  #
  #
  #
  # for(v in 1:numstudies){
  #
  #   study.names.list<-c(study.names.list,paste0("study",as.character(v)))
  #   betas.study.names.list<-c(betas.study.names.list,paste0("betas study ",as.character(v)))
  #   ses.study.names.list<-c(ses.study.names.list,paste0("ses study ",as.character(v)))
  # }
  #
  # # dimnames(betamatrix)<-list(dimnames(study.summary[[1]]$coefficients)[[1]], betas.study.names.list)
  # # dimnames(sematrix)<-list(dimnames(study.summary[[1]]$coefficients)[[1]], ses.study.names.list)
  # #
  # output.summary.text<-paste0("list(")
  #
  # for(u in 1:numstudies){
  #   output.summary.text<-paste0(output.summary.text,"study",as.character(u),"=study.summary[[",as.character(u),"]],"," ")
  # }
  #
  # output.summary.text.save<-output.summary.text
  # output.summary.text<-paste0(output.summary.text,"input.beta.matrix.for.SLMA=as.matrix(betamatrix),input.se.matrix.for.SLMA=as.matrix(sematrix))")
  #
  #
  # output.summary<-eval(parse(text=output.summary.text))
  #
  #
  # #######################END OF ANNOTATION CODE
  #
  # SLMA.pooled.ests.matrix<-matrix(NA,nrow<-numcoefficients,ncol=6)
  #
  #
  #
  #
  # if(!combine.with.metafor)
  # {
  #   return(output.summary)
  # }
  #
  # if(no.studies.valid)
  # {
  #   return(output.summary)
  # }
  #
  # #NOW ONLY WORKING WITH SITUATIONS WITH AT LEAST ONE VALID STUDY
  #
  # #IF combine.with.metafor == TRUE, FIRST CHECK THAT THE MODELS IN EACH STUDY MATCH
  # #IF THERE ARE DIFFERENT NUMBERS OF PARAMETERS THE ANALYST WILL
  # #HAVE TO USE THE RETURNED MATRICES FOR betas AND ses TO DETERMINE WHETHER
  # #COMBINATION ACROSS STUDIES IS POSSIBLE AND IF SO, WHICH PARAMETERS GO WITH WHICH
  # #ALSO DETERMINE WHICH STUDIES HAVE VALID DATA
  #
  # beta.matrix.for.SLMA<-as.matrix(betamatrix)
  # se.matrix.for.SLMA<-as.matrix(sematrix)
  #
  # #SELECT VALID COLUMNS ONLY (THERE WILL ALWAYS BE AT LEAST ONE)
  #
  # usecols<-NULL
  #
  # for(ut in 1:(dim(beta.matrix.for.SLMA)[2]))
  # {
  #   if(!is.na(beta.matrix.for.SLMA[1,ut])&&!is.null(beta.matrix.for.SLMA[1,ut]))
  #   {
  #     usecols<-c(usecols,ut)
  #   }
  # }
  #
  #
  # betamatrix.all<-beta.matrix.for.SLMA
  # sematrix.all<-se.matrix.for.SLMA
  #
  # betamatrix.valid<-beta.matrix.for.SLMA[,usecols]
  # sematrix.valid<-se.matrix.for.SLMA[,usecols]
  #
  # #CHECK FOR MATCHED PARAMETERS
  #
  # num.valid.studies<-as.numeric(dim(as.matrix(betamatrix.valid))[2])
  # coefficient.vectors.match<-TRUE
  #
  #
  #
  # if(num.valid.studies>1){
  #   for(j in 1:(num.valid.studies-1))
  #   {
  #     if(length(betamatrix.valid[,j])!=length(betamatrix.valid[,(j+1)]))coefficient.vectors.match<-FALSE
  #   }
  # }else{
  #   coefficient.vectors.match<-TRUE
  # }
  #
  #
  #
  #
  # if(!coefficient.vectors.match){
  #   cat("\n\nModels in different sources vary in structure\nplease match coefficients for meta-analysis individually\n")
  #   cat("nYou can use the DataSHIELD generated estimates and standard errors as the basis for a meta-analysis\nbut carry out the final pooling step independently of DataSHIELD using whatever meta-analysis package you wish\n\n")
  #   return(list(output.summary=output.summary))
  # }
  #
  #
  #
  # #IF combine.with.metafor == TRUE AND MODEL STRUCTURES MATCH ACROSS ALL STUDIES
  # #CREATE STUDY LEVEL META-ANALYSIS (SLMA) ESTIMATES FOR ALL PARAMETERS
  # #USING metafor() AND THREE APPROACHES TO SLMA: SLMA UNDER MAXIMUM LIKELIHOOD (SMLA-ML)
  # #SLMA UNDER RESTRICTED MAXIMUM LIKELIHOOD (SMLA-REML) AND USING FIXED EFFECTS (SLMA-FE)
  #
  # dimnames(SLMA.pooled.ests.matrix)<-list(dimnames(betamatrix.valid)[[1]],
  #                                         c("pooled.ML","se.ML","pooled.REML","se.REML","pooled.FE","se.FE"))
  #
  #
  #
  #
  #
  # for(p in 1:numcoefficients){
  #   rma.ML<-metafor::rma(yi=as.matrix(betamatrix.valid)[p,], sei=as.matrix(sematrix.valid)[p,], method="ML")
  #   rma.REML<-metafor::rma(yi=as.matrix(betamatrix.valid)[p,], sei=as.matrix(sematrix.valid)[p,], method="REML")
  #   rma.FE<-metafor::rma(yi=as.matrix(betamatrix.valid)[p,], sei=as.matrix(sematrix.valid)[p,], method="FE")
  #
  #   SLMA.pooled.ests.matrix[p,1]<-rma.ML$beta
  #   SLMA.pooled.ests.matrix[p,2]<-rma.ML$se
  #
  #   SLMA.pooled.ests.matrix[p,3]<-rma.REML$beta
  #   SLMA.pooled.ests.matrix[p,4]<-rma.REML$se
  #
  #   SLMA.pooled.ests.matrix[p,5]<-rma.FE$beta
  #   SLMA.pooled.ests.matrix[p,6]<-rma.FE$se
  #
  # }
  #
  #
  #
  # return(list(output.summary=output.summary, num.valid.studies=num.valid.studies,betamatrix.all=betamatrix.all,sematrix.all=sematrix.all, betamatrix.valid=betamatrix.valid,sematrix.valid=sematrix.valid,
  #             SLMA.pooled.ests.matrix=SLMA.pooled.ests.matrix))

}


# ds.coxSLMA
