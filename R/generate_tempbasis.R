#' Defining the temporal basis
#'
#' Use the data range and creates the basis functions of the temporal model
#'
#' @param data data.frame
#' @param tempBasis_fun string One of "re" (random effect) or "bs" (splines)
#' @param tempPeriod string Time unit (default to "%m" - month)
#' @param nSplines numeric Number of splines used to define the temporal basis functions
#' @param splinesType string One of "poly", "penalized", "cubic"
#' @export
generate_tempbasis <-   function(data, tempBasis_fun = "re", tempPeriod = "%m", nSplines = 3, splinesType = "poly"){

  SiteID = levels(data$ID)

  t = data$date

  BB = d2BB = iBB = NULL
  idxT = NULL
  j=0
  maxIdx = 0

  if(splinesType == "penalized"){
    for(i in tempPeriod){
      j = j +1
      if(j == 1){
        dat = data.frame(tTemp = as.numeric(as.Date(paste("2014-",format(t,"%m-%d"), sep = ""))))
        BBT <- mgcv::smoothCon(mgcv::s(tTemp,k=nSplines[j]),data=dat,knots=NULL)[[1]]
      }
      if(j == 2){
        dat = data.frame(tTemp = as.numeric(t))
        BBT <- mgcv::smoothCon(mgcv::s(tTemp,k=nSplines[j]),data=dat,knots=NULL, absorb.cons = F)[[1]]
      }
      if(j == 3){
        BBT = matrix(rep(1, length(t)), ncol=1)
      }

      BB = c(BB,list(BBT))
      idxt = 1:ncol(BBT$X) + maxIdx
      maxIdx = max(idxt)
      idxT = c(idxT,list(idxt))
    }
  }else{
    if(splinesType == "poly"){
      splineFun = function(...) splines::bs(intercept = F, ...)
      d2splineFun = function(...) splines2::dbs(intercept=F, derivs = 2, ...)
      isplineFun = function(...) splines2::ibs(intercept=F, ...)
    }
    if(splinesType == "cubic"){
      splineFun = function(...) splines::ns(intercept  =F, ...)
    }

    for(i in tempPeriod){
      j = j +1
      if(tempBasis_fun == "bs"){
        if(j == 1){
          if(tempPeriod[j] == "%m"){tTemp = as.Date(paste("2014-",format(t,"%m-%d"), sep = ""))}
          if(tempPeriod[j] == "%H"){tTemp = t}
          BBT = splineFun(tTemp, df = nSplines[j],
                          Boundary.knots = c(min(tTemp)-diff(range(tTemp))/(1*nSplines[j]),max(tTemp)+diff(range(tTemp))/(1*nSplines[j])))
          tTemp =as.numeric(tTemp)
          iBBT = apply(isplineFun(tTemp, df = nSplines[j],
                                  Boundary.knots = c(min(tTemp)-diff(range(tTemp))/(1*nSplines[j]),max(tTemp)+diff(range(tTemp))/(1*nSplines[j]))),2,max)
          d2BBT = colSums(d2splineFun(sort(unique(tTemp)), df = nSplines[j],
                                      Boundary.knots = c(min(tTemp)-diff(range(tTemp))/(1*nSplines[j]),max(tTemp)+diff(range(tTemp))/(1*nSplines[j]))))*diff(range(tTemp))
        }
        if(j == 2){
          BBT = splineFun(t, df = nSplines[j],
                          Boundary.knots = c(min(t)-diff(range(t))/(1*nSplines[j]),max(t)+diff(range(t))/(1*nSplines[j])))
          t = as.numeric(t)
          d2BBT = colSums(d2splineFun(sort(unique(t)), df = nSplines[j],
                                      Boundary.knots = c(min(t)-diff(range(t))/(1*nSplines[j]),max(t)+diff(range(t))/(1*nSplines[j]))))*diff(range(t))
          iBBT = apply(isplineFun(t, df = nSplines[j],
                                  Boundary.knots = c(min(t)-diff(range(t))/(1*nSplines[j]),max(t)+diff(range(t))/(1*nSplines[j]))),2,max)
        }
        if(j == 3){
          BBT = matrix(rep(1, length(t)), ncol=1)
        }
      }
      if(tempBasis_fun == "re"){
        BbaseT = as.numeric(unique(format(t, i)))
        BbaseTr = seq(min(BbaseT), max(BbaseT), by = 1)
        BBT = sapply(BbaseTr, function(x) as.numeric(as.numeric(format(t,i)) == x))
      }

      BB = c(BB,list(BBT))
      d2BB = c(d2BB,list(d2BBT))
      iBB = c(iBB,list(iBBT))
      idxt = 1:ncol(BBT) + maxIdx #(ncol(BB) - ncol(BBT))
      maxIdx = max(idxt)
      idxT = c(idxT,list(idxt))
    }
  }

  sptm_data = list(data = data, basis = BB, d2basis = d2BB, ibasis = iBB, list_idx = idxT, tempBasis_fun = tempBasis_fun, tempPeriod = tempPeriod, splinesType = splinesType)
  class(sptm_data) <- 'sptm'
  return(sptm_data)
}
