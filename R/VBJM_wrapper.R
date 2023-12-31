## HDJM wrapper functions ##

#' Simulated Longtidunal Data
#'
#' This dataset contains longitudinal outcomes.
#'
#' \itemize{
#'   \item ID patient ID
#'   \item item types of longitudinal outcome
#'   \item years measurement timepoints
#'   \item value measurements
#' }
#'
#' @name LongData
#' @docType data
#' @author Jiehuan Sun \email{jiehuan.sun@gmail.com}
#' @keywords data
#' @usage data(VBJMdata)
#' @format A data frame with 48700 rows and 4 variables
NULL

#' Simulated Survival Data
#'
#' This dataset contains survival outcome.
#'
#' \itemize{
#'   \item ID patient ID
#'   \item fstat censoring indicator
#'   \item ftime survival time
#'   \item x baseline covariates
#' }
#' @name SurvData
#' @docType data
#' @author Jiehuan Sun \email{jiehuan.sun@gmail.com}
#' @keywords data
#' @usage data(VBJMdata)
#' @format A data frame with 100 rows and 4 variables
NULL


#' control_list
#'
#' This list contains a list of parameters specifying the joint model.
#'
#' \itemize{
#'   \item ID_name the variable name for the patient ID in both
#'   longitudinal data and survival data.
#'   \item item_name the variable name for the longitudinal outcomes
#'   in the longitudinal data.
#'   \item value_name the variable name for the longitudinal measurements
#'   in the longitudinal data.
#'   \item time_name the variable name for the measurement timepoints in the
#'   longitudinal data.
#'   \item fix_cov a set of variables names indicating the covariates of
#'   fixed-effects in the longitudinal submodel.
#'   If NULL, not baseline covariates are included.
#'   \item random_cov a set of variables names indicating the covariates of
#'   random-effects in the longitudinal submodel.
#'   If NULL, not baseline covariates are included.
#'   \item FUN a function specifying the time-related basis functions in
#'   the longitudinal submodel.
#'   \item ran_time_ind a vector of integers specifying which
#'   time-related basis functions are also included with random-effects in
#'   the longitudinal submodel.
#'   \item surv_time_name the variable name for the survival time
#'   in the survival data.
#'   \item surv_status_name the variable name for the censoring indicator
#'   in the survival data.
#'   \item surv_cov a set of variables names specifying the baseline covariates
#'   in the survival submodel.
#'   \item n_points an integer indicating the numebr of nodes being used in
#'   the Gaussian quadrature.
#' }
#' @name control_list
#' @docType data
#' @author Jiehuan Sun \email{jiehuan.sun@gmail.com}
#' @keywords data
NULL




#' @noRd
prep_data <- function(LongData=NULL, SurvData = NULL,
                      control_list=NULL, marker.name=NULL){

    ### control_list
    ID_name = control_list$ID_name
    item_name = control_list$item_name
    value_name = control_list$value_name
    time_name = control_list$time_name
    fix_cov = control_list$fix_cov
    random_cov = control_list$random_cov
    FUN = control_list$FUN
    ran_time_ind=control_list$ran_time_ind
    surv_time_name = control_list$surv_time_name
    surv_status_name = control_list$surv_status_name
    surv_cov =  control_list$surv_cov
    n_points =  control_list$n_points

    ###
    data.list = list()
    para.list = list()

    flex_time = FUN(1)
    fix_est_name = c("intercept", fix_cov, colnames(flex_time))
    rand_est_name = c("intercept",random_cov, colnames(flex_time)[ran_time_ind])
    surv_est_name = surv_cov

    ## run LME to initiate the parameters in Longitudinal submodel
    ## Y_i
    uni_ID = SurvData[,ID_name]
    if(is.null(marker.name)){
        marker.name = unique(LongData[,item_name])
    }

    Y.list = lapply(uni_ID, function(i){
        data.tmp = LongData[LongData[,ID_name]==i,]
        lapply(marker.name,function(x){
            matrix(data.tmp[data.tmp[,item_name]==x,value_name],ncol=1)
            # matrix(data.tmp$value[data.tmp$item==x],ncol=1)
        })
    })

    Y.list = do.call(rbind, Y.list)

    ## X_i, Z_i
    X.list = lapply(uni_ID, function(i){
        data.tmp =  LongData[LongData[,ID_name]==i,]

        lapply(marker.name,function(x){
            if(!is.null(fix_cov)){
                time_mat.tmp = FUN(data.tmp[data.tmp[,item_name]==x,time_name])
                cov.tmp = SurvData[SurvData[,ID_name]==i, fix_cov, drop=FALSE]
                cov.tmp = as.matrix(cov.tmp[rep(1,nrow(time_mat.tmp)),])
                as.matrix(cbind(1, cov.tmp, time_mat.tmp ))
            }else{
                as.matrix(cbind(1, FUN(data.tmp[data.tmp[,item_name]==x,time_name])))
            }

            #cbind(1,data.tmp$years[data.tmp$item==x])
        })
    })
    X.list = do.call(rbind, X.list)


    Z.list = lapply(uni_ID, function(i){
        data.tmp =  LongData[LongData[,ID_name]==i,]

        lapply(marker.name,function(x){
            if(!is.null(random_cov)){
                time_mat.tmp = FUN(data.tmp[data.tmp[,item_name]==x,time_name])
                cov.tmp = SurvData[SurvData[,ID_name]==i, random_cov, drop=FALSE]
                cov.tmp = as.matrix(cov.tmp[rep(1,nrow(time_mat.tmp)),])
                as.matrix(cbind(1, cov.tmp, time_mat.tmp[,ran_time_ind,drop=FALSE]))
            }else{
                as.matrix(cbind(1,  FUN(data.tmp[data.tmp[,item_name]==x,time_name])[,ran_time_ind,drop=FALSE]))
            }

            #cbind(1,data.tmp$years[data.tmp$item==x])
        })
    })
    Z.list = do.call(rbind, Z.list)


    ## X_i(T_i),  z_i(T_i)
    X_T.list = lapply(uni_ID, function(i){
        T_i = SurvData[SurvData[,ID_name] == i,surv_time_name]
        #T_i = SurvData$ftime[SurvData$ID==i]
        data.tmp =  LongData[LongData[,ID_name]==i,]

        lapply(marker.name,function(x){
            vv = SurvData[SurvData[,ID_name]==i, fix_cov, drop=FALSE]
            matrix(c(1, as.numeric(vv[1,]), as.numeric(FUN(T_i))), ncol=1)
            #matrix(c(1,T_i),ncol=1)
        })
    })
    X_T.list = do.call(rbind, X_T.list)

    Z_T.list = lapply(uni_ID, function(i){
        T_i = SurvData[SurvData[,ID_name] == i,surv_time_name]
        #T_i = SurvData$ftime[SurvData$ID==i]
        data.tmp =  LongData[LongData[,ID_name]==i,]
        lapply(marker.name,function(x){
            vv = SurvData[SurvData[,ID_name]==i, random_cov, drop=FALSE]
            matrix(c(1,as.numeric(vv[1,]) ,as.numeric(FUN(T_i))[ran_time_ind]), ncol=1)
            #matrix(c(1,T_i),ncol=1)
        })
    })
    Z_T.list = do.call(rbind, Z_T.list)

    ## X_i(t) , z_i(t)
    ## this depends on the number of legendre Gaussian quadrature points
    Gauss.point  = gauss.quad(n_points)
    # \int_0^{T_i} f(t)dt
    # t_node = Gauss.point$nodes *(Ti/2) + Ti/2
    # w_node = Gauss.point$weights
    # Ti/2 * sum(w_node * f(t_node))

    X_t.list = lapply(uni_ID, function(i){
        Ti = SurvData[SurvData[,ID_name] == i,surv_time_name]
        t_node = Gauss.point$nodes *(Ti/2) + Ti/2
        data.tmp =  LongData[LongData[,ID_name]==i,]
        time_mat.tmp = FUN(t_node)

        lapply(marker.name,function(x){
            if(!is.null(fix_cov)){
                cov.tmp = SurvData[SurvData[,ID_name]==i, fix_cov, drop=FALSE]
                cov.tmp = as.matrix(cov.tmp[rep(1,nrow(time_mat.tmp)),])
                as.matrix(cbind(1, cov.tmp, time_mat.tmp ))
            }else{
                as.matrix(cbind(1,  time_mat.tmp ))
            }
        })
    })

    X_t.list = do.call(rbind, X_t.list)

    Z_t.list = lapply(uni_ID, function(i){
        Ti = SurvData[SurvData[,ID_name] == i,surv_time_name]
        t_node = Gauss.point$nodes *(Ti/2) + Ti/2
        data.tmp =  LongData[LongData[,ID_name]==i,]
        time_mat.tmp = FUN(t_node)

        lapply(marker.name,function(x){

            if(!is.null(random_cov)){
                cov.tmp = SurvData[SurvData[,ID_name]==i, random_cov, drop=FALSE]
                cov.tmp = as.matrix(cov.tmp[rep(1,nrow(time_mat.tmp)),])
                as.matrix(cbind(1, cov.tmp, time_mat.tmp[,ran_time_ind,drop=FALSE] ))
            }else{
                as.matrix(cbind(1,  time_mat.tmp[,ran_time_ind,drop=FALSE] ))
            }
        })
    })

    Z_t.list = do.call(rbind, Z_t.list)


    w_node.list = lapply(uni_ID, function(i){
        Ti = SurvData[SurvData[,ID_name] == i,surv_time_name]
        Gauss.point$weights*Ti/2
    })
    t_node.list = lapply(uni_ID, function(i){
        Ti = SurvData[SurvData[,ID_name] == i,surv_time_name]
        Gauss.point$nodes *(Ti/2) + Ti/2
    })

    ## covariates in survival submodel
    W = as.matrix(SurvData[,surv_cov,drop=FALSE])

    data.list[["Y"]] = Y.list
    data.list[["X"]] = X.list
    data.list[["X_t"]] = X_t.list
    data.list[["X_T"]] = X_T.list
    data.list[["Z"]] = Z.list
    data.list[["Z_t"]] = Z_t.list
    data.list[["Z_T"]] = Z_T.list
    data.list[["W"]] = W
    data.list[["GQ_w"]] = w_node.list
    data.list[["GQ_t"]] = t_node.list
    data.list[["ftime"]] = SurvData[,surv_time_name]
    data.list[["fstat"]] = SurvData[,surv_status_name]

    list(data.list=data.list, uni_ID=uni_ID, marker.name=marker.name,
         fix_est_name=fix_est_name,rand_est_name=rand_est_name,
         surv_est_name=surv_est_name)
}

#' @noRd
VBJM_init <- function(LongData=NULL, SurvData = NULL, control_list=NULL,
                      marker.name=NULL){

    ###
    data.list = prep_data(LongData=LongData, SurvData = SurvData,
                          control_list=control_list,marker.name=marker.name)
    #uni_ID = data.list$uni_ID
    marker.name = data.list$marker.name
    fix_est_name = data.list$fix_est_name
    rand_est_name = data.list$rand_est_name
    surv_est_name = data.list$surv_est_name
    data.list = data.list$data.list

    ## run LME to initiate the parameters in Longitudinal submodel
    para.list = list()
    beta.list = list()
    mu.list = list()
    V.list = list()
    Sigma.list = list()
    sig.vec = rep(NA, length(marker.name))
    alpha.vec = rep(NA,length(marker.name))

    for(i in seq_along(marker.name)){
        # i = 1
        # print(i)
        fitLME = init_LME(data.list$Y[,i], data.list$X[,i],
                          data.list$Z[,i], 100, 1e-4)
        beta.list[[i]] = fitLME$beta
        mu.list[[i]] = fitLME$mu
        sig.vec[i] = fitLME$sig2
        Sigma.list[[i]] = fitLME$Sigma
        V.list[[i]] = fitLME$V
        alpha.vec[i] = 0
    }

    Sigma = as.matrix(bdiag(Sigma.list))

    V.list = do.call(cbind, V.list)
    V.list = lapply(1:nrow(V.list), function(i){
        as.matrix(bdiag(  V.list[i,] ))
    })

    mu.list = do.call(cbind, mu.list)

    ## initiate the parameters in Survival submodel

    fitSURV = survreg(Surv(data.list[["ftime"]],data.list[["fstat"]] ) ~ data.list[["W"]])

    theta = exp(fitSURV$coefficients[1])
    lambda = 1/fitSURV$scale
    gamma = -fitSURV$coefficients[2:(1+length(surv_est_name))]/fitSURV$scale

    ###
    para.list[["mu"]] = mu.list
    para.list[["V"]] = V.list
    para.list[["Sigma"]] = Sigma
    para.list[["sig2"]] = sig.vec

    para.list[["beta"]] = beta.list
    para.list[["weib"]] = c(lambda, theta)
    para.list[["gamma"]] = gamma
    para.list[["alpha"]] = alpha.vec

    list(data.list=data.list, para.list=para.list,
         marker.name=marker.name, fix_est_name=fix_est_name,
         rand_est_name=rand_est_name, surv_est_name=surv_est_name
    )
}


#' @noRd
VBJM_get_summary <- function(init_list=NULL, res=NULL){

    marker.name = init_list$marker.name
    fix_est_name = init_list$fix_est_name
    #rand_est_name = init_list$rand_est_name
    surv_est_name = init_list$surv_est_name

    beta.list = lapply(seq_along(marker.name), function(i){
        beta = as.numeric(res$beta[[i]])
        coef_name = paste(marker.name[i],"_fix_",fix_est_name,sep="")
        names(beta) = coef_name
        beta
    })

    gamma = as.numeric(res$gamma)
    names(gamma) = paste("Surv_gamma_",surv_est_name,sep="")

    alpha = as.numeric(res$alpha)
    names(alpha) = paste(marker.name,"_alpha", sep="")

    weib = as.numeric(res$weib)
    names(weib) = c("Weibull_shape","Weibull_scale")

    para =c(do.call(c, beta.list),  gamma, alpha, weib)

    #res$H = -res$H
    #diag(res$H) = diag(res$H) + 1e-6
    #cov = solve(res$H)
    cov = -pinv(res$H)
    se = round(sqrt(diag(cov)),4)[1:length(para)]

    res_summary = data.frame(Estimate=para, SE=se,
                             para-1.96*se, para+1.96*se)

    se_weib = se[c(length(se)-1, length(se))]
    se_log_weib = sqrt(se_weib^2 / weib^2)

    ci_weib_1 = exp(log(weib[1]) + c(-1.96, 1.96) *se_log_weib[1])
    ci_weib_2 = exp(log(weib[2]) + c(-1.96, 1.96) *se_log_weib[2])
    res_summary[c(length(se)-1, length(se)),3:4] = rbind(ci_weib_1, ci_weib_2)

    colnames(res_summary)[3:4] = c("95%CI_lower","95%CI_upper")
    res_summary
}


#' The function to fit VBJM.
#'
#' The function is used to fit joint models using variational inference algorithm.
#'
#' @param LongData a data frame containing the longitudinal data
#' (see \code{\link{LongData}}).
#' @param SurvData a data frame containing the survival data
#' (see \code{\link{SurvData}}).
#' @param marker.name a vector indicating which set of longitudinal biomarkers
#' to be analyzed. If NULL, all biomarkers in LongData will be used.
#' @param control_list a list of parameters specifying the joint model
#' (see \code{\link{control_list}}).
#' @param maxiter the maximum number of iterations.
#' @param eps threshold for convergence.
#'
#' @return return a data frame with estimates, standard errors, and
#' 95\% CIs for each of the following parameters, where VAR indicates
#' the corresponding variable name.
#' \item{VAR_alpha}{the parameters for the effects of biomarkers
#' in the survival submodel, where VAR indicates the names
#' for the biomarkers.}
#' \item{Weibull_shape}{the shape parameter in the
#'  Weibull baseline hazard in the survival submodel.}
#' \item{Weibull_scale}{the scale parameter in the
#'  Weibull baseline hazard in the survival submodel.}
#'  \item{Surv_gamma_VAR}{the parameters for the effects of baseline covariates
#' in the survival submodel.}
#' \item{VAR_fix}{the parameters for the fixed-effects
#' in the longitudinal submodel.}
#'
#' @references Jieqi Tu and Jiehuan Sun (2023). "Gaussian variational
#' approximate inference for joint models of longitudinal biomarkers
#'  and a survival outcome". Statistics in Medicine, 42(3), 316-330.
#'
#' @examples
#' data(VBJMdata)
#' flex_time_fun <- function(x=NULL){
#'     xx = matrix(x, ncol = 1)
#'     colnames(xx) = c("year_l")
#'     xx
#' }
#' ran_time_ind = 1 ## random time-trend effects
#' control_list = list(
#'   ID_name = "ID", item_name = "item",
#'   value_name = "value",  time_name = "years",
#'   fix_cov = NULL, random_cov = NULL,
#'   FUN = flex_time_fun, ran_time_ind=ran_time_ind,
#'   surv_time_name = "ftime",  surv_status_name = "fstat",
#'   surv_cov = "x", n_points = 5
#' )
#' \donttest{
#' ## takes about one minute.
#' res = VBJM_fit(LongData=LongData, SurvData=SurvData,
#'                control_list=control_list)
#' }
#'
VBJM_fit <- function(LongData = NULL, SurvData = NULL,marker.name = NULL,
                     control_list = NULL,  maxiter=100, eps=1e-4){

    ## intiate the algorithm by re-formatting the dataset
    message("initialization ...")
    init_list = VBJM_init(LongData=LongData, SurvData = SurvData,
                          control_list=control_list, marker.name=marker.name)
    ## running lasso
    message("running VBJM ...")
    fitVBJM = VBJM(init_list$data.list,  init_list$para.list,
                   maxiter=maxiter, eps=eps)
    res_VBJM = VBJM_get_summary(init_list=init_list, res=fitVBJM)
    res_VBJM
}






