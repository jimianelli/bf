
.THEME    <<- ggthemes::theme_few(base_size = 12, base_family = "")
plot_fmort <- function(M,runs=1:5,wrap=FALSE){
mdf <- data.frame(NULL)
for (i in runs) {
	A <- M[[i]]
  df <- data.frame(Fmort=A$Fmort)
  df$model = names(M)[i]
  df$Year= 1978:(1977+length(A$SSB))
  df$lb = A$Fmort-2*A$Fmort.sd
  df$ub = A$Fmort+2*A$Fmort.sd
  mdf <- rbind(mdf,df)
}
if (wrap)
  p <- ggplot(mdf,aes(x=Year,y=Fmort,fill=model)) + .THEME + geom_ribbon(aes(ymin=lb,ymax=ub),alpha=.2) + geom_line(aes(col=model)) + facet_grid(.~model)
else
  p <- ggplot(mdf,aes(x=Year,y=Fmort,fill=model)) + .THEME + geom_ribbon(aes(ymin=lb,ymax=ub),alpha=.2) + geom_line(aes(col=model))
p <- p + expand_limits(y = 0)
return(p)
}

plot_ssb <- function(M,runs=1:5,wrap=FALSE){
mdf <- data.frame(NULL)
for (i in runs) {
	A <- M[[i]]
  df <- data.frame(SSB=A$SSB)
  df$model = names(M)[i]
  df$Year= 1978:(1977+length(A$SSB))
  df$lb = A$SSB-2*A$SSB.sd
  df$ub = A$SSB+2*A$SSB.sd
  mdf <- rbind(mdf,df)
}
if (wrap)
  p <- ggplot(mdf,aes(x=Year,y=SSB,fill=model)) + .THEME + geom_ribbon(aes(ymin=lb,ymax=ub),alpha=.2) + geom_line(aes(col=model)) + facet_grid(model~.)
else
  p <- ggplot(mdf,aes(x=Year,y=SSB,fill=model)) + .THEME + geom_ribbon(aes(ymin=lb,ymax=ub),alpha=.2) + geom_line(aes(col=model))
p <- p + expand_limits(y = 0)
return(p)
}

plot_rec <- function(M,runs=1:5,wrap=FALSE){
mdf <- data.frame(NULL)
for (i in runs) {
	A <- M[[i]]
  df <- data.frame(rec=A$rec)
  df$model = names(M)[i]
  df$Year= 1978:(1977+length(A$rec))
  df$lb=(A$rec/exp(2.*sqrt(log(1+(A$rec.sd^2)/(A$rec^2)))));
  df$ub=(A$rec*exp(2.*sqrt(log(1+(A$rec.sd^2)/(A$rec^2)))));
  #df$lb = A$rec-2*A$rec.sd
  #df$ub = A$rec+2*A$rec.sd
  mdf <- rbind(mdf,df)
}
if (wrap)
  p <- ggplot(mdf,aes(x=Year,y=rec,fill=model)) + .THEME + geom_errorbar(aes(ymin=lb,ymax=ub,color=model),position="dodge2",alpha=.8) + geom_bar(aes(fill=model),position="dodge",stat="identity",alpha=.5) + facet_grid(.~model)
else
  p <- ggplot(mdf,aes(x=Year,y=rec,fill=model)) + .THEME + geom_errorbar(aes(ymin=lb,ymax=ub,color=model),position="dodge2",alpha=.8) + geom_bar(aes(fill=model),position="dodge2",stat="identity",alpha=.5)
p <- p + ylab("Recruitment")
return(p)
}

plot_ind<- function(M,runs=c(2:5),wrap=FALSE){
mdf <- data.frame(NULL)
for (i in runs) {
	A <- M[[i]]
  df <- data.frame(pred=A$pred_trend)
  df$obs =A$obs_trend
  df$model = names(M)[i]
  df$Year= A$yrs_trend
  #df$lb = A$rec-2*A$rec.sd
  #df$ub = A$rec+2*A$rec.sd
  mdf <- rbind(mdf,df)
}
  p <- ggplot(mdf,aes(x=Year,y=obs)) + .THEME +  
    geom_bar(position="dodge2",stat="identity") + facet_grid(model~.) + 
    geom_line(data=mdf,aes(x=Year,y=pred,color=model),size=2)
return(p)
}


tab_fit <- function(M,mod_scen=NULL){
  df <- NULL
  if (is.null(mod_scen)) mod_scen=1:length(M)
  for (ii in mod_scen) {
    x         <- M[[ii]]
    cat_like  <- round(x$like[1],2); names(cat_like) <- paste0("Catch NLL")
    fac_like  <- round(x$like[2],2); names(fac_like) <- paste0("Fish Age NLL")
    ind_like  <- round(x$like[3],2); names(ind_like) <- paste0("Index NLL")
    r_like    <- round(x$like[4],2); names(r_like) <- paste0("Recruit penalty NLL")
    q_like    <- round(x$like[5],2); names(q_like) <- paste0("q penalty NLL")
    f_prior  <- round(sum(x$prior[1]),2); names(f_prior) <- paste0("F Prior")
    m_prior  <- round(sum(x$prior[2]),2); names(m_prior) <- paste0("M Prior")
    s_prior  <- round(sum(x$prior[3]),2); names(s_prior) <- paste0("selectivity Prior")
    q_prior  <- round(sum(x$prior[4]),2); names(q_prior) <- paste0("q Prior")
    rw_prior  <- round(sum(x$prior[5]),2); names(rw_prior) <- paste0("rand-walk Prior")
    Priors   <- sum(x$prior); names(Priors) <- "Priors"
    Data   <- sum(x$like) ; names(Data) <- "Data and penalties"
    Total <- Data+Priors ; names(Total) <- "Total"
    v      <- c(  cat_like,ind_like, r_like, q_like, 
                 fac_like, f_prior, m_prior, s_prior, q_prior, rw_prior, Priors,
                 Data, Total)
    df     <- cbind(df, v)
  }
  df <- data.frame(rownames(df), df, row.names = NULL)
  names(df) <- c("Component",names(M[mod_scen]))
  return(df)
}
