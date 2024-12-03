#' Create spatio temporal data
#'
#' Create fake data following a cluster of spatio-temporal signals
#'
#' @param min_date date minimum data for the simulated data
#' @param max_date date maximum data for the simulated data
#' @param frequency string Observation frequency. Default to "month".
#' @param nSite numeric Number of locations.
#' @param nLandUse numeric Number of different landuses.
#' @param nBasis numeric Number of temporal basis functions.
#' @param betaTrue numeric matrix Range of variation. Default to -2 to 2.
#' @param parameters list Uncertainty parameters. Default to sigma (1), tau (0.5), phi (1).
#' @param missingMeasures list Frequency and type of missing data. Default to no missing measures.
#' @export
generate_sptmixtdata <-  function(min_date = as.Date("01-01-2009", format = "%d-%m-%Y"),
           max_date = as.Date("31-12-2010", format = "%d-%m-%Y"),
           frequency = "month",
           nSite = 8,
           nLandUse = 2,
           nBasis = 2,
           betaTrue = NULL,
           parameters = NULL,
           missingMeasures = list(status = F, ratio = 0, type="MNAR")
  ){

    date = seq(min_date,max_date,by = frequency)
    seasonality = array(NA, dim = c(length(date),2,nLandUse))
    shift = sample(-6:6, size = nLandUse, replace = F)
    if(frequency == "month"){
      for(i in 1:nLandUse){
        seasonality[,,i] = cbind(cos(((as.numeric(format(date, "%m")) + shift[i])*2*3.146 )/ 12),sin(((as.numeric(format(date, "%m")) + shift[i]))*2*3.146 / 12)^2)
      }
    }
    if(frequency =="week"){
      for(i in 1:nLandUse){
        seasonality[,,i] = cbind(cos(((as.numeric(format(date, "%W")) + shift[i])*2*3.146) / 52),sin(((as.numeric(format(date, "%W")) + shift[i])*2*3.146) / 52)^2)
      }
    }

    if(is.null(betaTrue)){
      betaTrue = matrix(runif((nBasis+1)*nLandUse,-2,2),ncol=nBasis+1)
    }
    if(is.null(parameters)){
      parameters = list(sigma = 0.1, tau = 0.5, phi = 1)
    }

    locations = matrix(20*runif(nSite*2),ncol=2)
    LU.comp = dirmult::rdirichlet(nSite,alpha = rep(0.5,nLandUse))
    LU.true = t(apply(LU.comp,1,function(x) rmultinom(1,1,x)))

    sim.raw.cov = data.frame(ID = as.factor(1:nSite), Longitude = locations[,1], Latitude = locations[,2])
    LU.df = as.data.frame(LU.comp)
    names(LU.df) <- paste0("LU",1:nLandUse)
    LU.df$ID = 1:nSite

    sim.raw.cov = merge(sim.raw.cov, LU.df, by = "ID")
    sim.raw.obs = data.frame(date = rep(date, each= nSite), ID = as.factor(rep(1:nSite, length(date))))
    sim.raw.obs$obs = NA

    X.true = LU.true
    D = as.matrix(dist(sim.raw.cov[,c("Longitude","Latitude")]))
    for(i in 1:nSite){
      sim.raw.obs$obs[sim.raw.obs$ID == i] = betaTrue[which.max(X.true[i,]),1] +
        seasonality[,1:nBasis, which.max(X.true[i,])] %*% matrix(betaTrue[which.max(X.true[i,]),2:(2+nBasis-1)], ncol=1)
    }

    SpatialRF = mvtnorm::rmvnorm(1, mean=rep(0,nSite), sigma = parameters$tau*exp(-parameters$phi*D))

    for(j in date){
      sim.raw.obs$obs[sim.raw.obs$date == j] = sim.raw.obs$obs[sim.raw.obs$date == j] + SpatialRF + rnorm(length(diag(D)),0,parameters$sigma)
    }
    if(missingMeasures$status){
      if(missingMeasures$type == "MNAR"){
        keepDates = sort(sample(x = 1:length(date), size = round(length(date)*(1-missingMeasures$ratio))))
        sim.raw.obs = sim.raw.obs[sim.raw.obs$date %in% date[keepDates],]}
      if(missingMeasures$type == "MAR"){
        keepObs = sort(sample(x = 1:nrow(sim.raw.obs), size = round(nrow(sim.raw.obs)*(1-missingMeasures$ratio))))
        sim.raw.obs = sim.raw.obs[keepObs,]}
    }
    return(list(data = sim.raw.obs, X = sim.raw.cov, coordinates = locations, luComp = LU.df,
                betaTrue = betaTrue, grTrue = LU.true, parTrue = parameters))
  }
