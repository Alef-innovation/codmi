
####################################################################
##########           Functions definition               ############
{
  #CoDMI parameters: chosen values and default settings
  {
    epsilon=1            # tolerance (in time units)
    maxiter=1000         # maximum number of allowed iterations
    convcrit=0           # convergence criterion. 0: maximum absolute error
    init_e=3             # life expectancies inits. 3: computed on DoD+Cen cases
    compfinpoint=TRUE    # survival distribution: completed
    conf_lev=0.95        # confidence level for extended Greenwood's formula
    cenadj=FALSE         # final adjustment for censoring
    type.alpha=4         # prob(DoD/Cen). 4: computed through hazard rate funtions  
    exo_alpha=NULL       # prob(DoD/Cen): exogenously inputed
  }
  
  #############    CoDMI  algorithm     ##############
  codmi=function(td_dodcen, t_cov,
                 init_e=3, epsilon=1, maxiter=1000, 
                 cenadj=FALSE,exo_alpha=NULL,
                 conf_lev=0.95,...){
    theta=t_cov
    df_t=td_dodcen
    nt=length(t_cov)
    
    ###  Initialization step    #########
    #initialize life expectancies
    if(init_e==1) e_old=0
    if(init_e==2) e_old=theta     
    if(init_e==3){
      pdfz=KM_AV(df_t,compfinpoint=compfinpoint)
      e_old=sapply(1:nt,function(x) Exp_t_theta(pdfz,theta[x]))-theta
    }
    if(init_e==4) e_old=max(c(df_t$t,theta))*runif(nt)-theta
    # initialize status indicators
    d_j=rep(1,nt)
    
    ### Iterative loop   #########
    out=matrix(NA,1,5)
    iter=0
    mean_err=1000
    while ( (mean_err>epsilon) & (iter < maxiter) ){
      iter=iter+1
      
      # new complete data set
      data.all=rbind(df_t,
                     data.frame(t=theta+e_old,d=d_j))
      data.all=data.all[order(data.all$t),]
      
      ### Estimation step (survival function estimated)  #########
      pdf.all=KM_AV(data.all,compfinpoint=compfinpoint)
      ll=pdf.all$ll_sum[1]
      ### Expectation step (expected lifetimes computed)  #########
      tau = sapply(1:nt,function(x) Exp_t_theta(pdf.all,theta[x]))
      
      e_hat=tau-theta
      
      mean_e_hat=mean(e_hat)
      mean_e_old=mean(e_old)
      
      if(convcrit==0){mean_err=max(abs(e_hat-e_old))}
      if(convcrit==1){mean_err=mean(abs(e_hat-e_old))}
      if(convcrit==2){mean_err=abs(mean_e_hat-mean_e_old)}
      
      out=rbind(out,
                cbind(iter,mean_e_old,mean_e_hat,mean_err,loglik=ll))
      e_old=e_hat
    }
    mean_iter=out[2:nrow(out), ]
    mean_fin=data.frame(t(out[nrow(out), ]))
    
    ### Final censoring adjustment  #########
    theta_cen=NA
    if(cenadj==TRUE){
      e_hat1=e_hat
      data.all.rev=data.all
      data.all.rev$d=1-data.all.rev$d
      pdf.all0=KM_AV(data.all.rev,compfinpoint=compfinpoint)
      
      tau0 = sapply(1:nt,function(x) Exp_t_theta(pdf.all0,theta[x]))
      e_hat0=tau0-theta
      
      # compute alpha_j
      alpha_j=compute_alpha(df_t,theta,type.alpha)
      
      # case exo_alpha!=NULL
      if(!is.null(exo_alpha)){
        if (length(exo_alpha)!=nt) 
          stop("exo_alpha has a different length than theta")
        alpha_j=exo_alpha
      }
      
      d_j[alpha_j<0.5]=0
      e_hat=e_hat1
      e_hat[alpha_j<0.5]=e_hat0[alpha_j<0.5]
      theta_cen=data.frame(theta,e_hat,e_hat1,e_hat0,alpha_j)
      mean_fin$mean_e_hat=mean(e_hat)
    }
    
    
    table_fin <- data.frame(j=1:nt,
                            theta=theta,
                            e_hat=e_hat,
                            tau_hat=theta+e_hat,
                            delta=d_j)
    niter <- iter
    if (niter >= maxiter){
      converge = FALSE
      convmsg="Algorithm failed to converge"
    }else{
      converge = TRUE 
      convmsg="Convergence criterion met"
    }
    print(convmsg)
    
    ###   Adjusted Variance (extended Greenwood's formula)  #########
    
    # CoDMI output (w_hat)
    data.all=rbind(df_t,
                   data.frame(t=table_fin$tau_hat,
                              d=table_fin$delta))
    data.all=data.all[order(data.all$t),]
    
    # Compute q_i probabilities  
    df=data.all
    pdf=KM_AV(df,adjvar=FALSE,delcen=FALSE,delties=FALSE)
    
    # Reverse KM 
    df_rev=df
    df_rev$d=1-df_rev$d
    pdf_rev=KM_AV(df_rev,adjvar=FALSE,delcen=FALSE,delties=FALSE)
    
    # Compute Q_i probabilities 
    
    q_ij=compute_q_ij(pdf,theta)
    q_ij_rev=compute_q_ij(pdf_rev,theta)
    
    d_ij=matrix(table_fin$delta,nrow(q_ij),ncol(q_ij),byrow = TRUE)
    q_ij_w=d_ij*q_ij+(1-d_ij)*q_ij_rev
    Q_i=rowSums(q_ij_w)
    Q2_i=rowSums(q_ij_w**2)
    
    pdf$Q=Q_i  
    pdf$Q2=Q2_i
    
    # Compute Deltas and reorder  
    
    pdf$Delta=0
    pdf[pdf$t %in% table_fin$tau_hat,]$Delta=1
    pdf=pdf[order(pdf$t),]
    
    pdfAV=KM_AV(pdf,adjvar = TRUE)
    pdfC =KM_AV(pdf)
    
    # extended c.i. 
    zalpha <- qnorm(1-(1-conf_lev)/2)
    coeff=pmax(pdfC$S,pdfAV$S)/pdfAV$S
    ci_lower_AV=exp(log(pdfC$S)-zalpha*pdfAV$cv*coeff)
    ci_upper_AV=exp(log(pdfC$S)+zalpha*pdfAV$cv*coeff)
    # restrict
    ci_lower_AV[ci_lower_AV < 0]=0
    ci_upper_AV[ci_upper_AV > 1]=1
    
    # classical c.i. (Greenwood)
    ci=comp_ciKM(pdfC)
    
    pdf_fin=data.frame(pdfAV,
                       ci_lower=ci$ci_lower_S,
                       ci_upper=ci$ci_upper_S,
                       ci_lower_AV,
                       ci_upper_AV)
    # list of outputs
    l_out=list(mean_iter=mean_iter,
               table_fin=table_fin,
               converge=converge,
               niter=niter,
               theta_cen=theta_cen,
               mean_fin=mean_fin,
               pdf_fin=pdf_fin)
    return(l_out)
  }
  
  #############     Exp_t_theta function     ##############
  # Compute expected lifetime given t>theta
  # (if theta>max(t) then Exp_t_theta=max(t))
  # input: df: data.frame(t,p), theta: scalar
  # output: Exp_t_theta (scalar)
  Exp_t_theta=function(df, theta){
    c=(df$t>theta)
    out=ifelse(sum(c)>0,
               sum(df$t[c]*df$p[c])/sum(df$p[c]),
               max(df$t))
    return(out)
  }
  
  ###############   KM_AV  function     #########################
  # Compute KM estimates with possible adjustment for variance
  # input: df: data.frame(t,d, and Delta, Q, Q2 if adjvar=T)
  #        adjvar=T/F        T: variance adjusted for DoC cases
  #        delcen=T/F        T: censored values not reported
  #        delties=T/F       T: last tied value taken
  #        compfinpoint=T/F  T: incomplete distribution completed
  # output: data.frame(df,R,Y,nCen,S,sd_S,cv,ll)
  KM_AV=function(df,adjvar=FALSE,
                 delcen=TRUE,delties=TRUE,compfinpoint=TRUE){
    
    #completeness controlled
    if (adjvar) {
      colsel=c("t","d","Delta","Q","Q2")
      if (sum(colnames(df) %in% colsel)!=length(colsel)) 
        stop("data.frame doesn't contain t, d, Delta, Q, Q2")
      df=df[,colsel]  # input selected
    }else{
      colsel=c("t","d")
      if (sum(colnames(df) %in% colsel)!=length(colsel))
        stop("data.frame doesn't contain t, d")
      df=df[,colsel]  # input selected
    }
    # ordering
    df=df[order(df$t,-df$d),]  #in case of tie, Cen after DoD 
    
    zn=1e-8
    n=nrow(df)
    R=rep(NA,n)     # R, number of subjects at risk 
    S=rep(NA,n)     # S_hat_KM
    cv2=rep(NA,n)   # cv2_S_hat_KM
    ties=rep(0,n)    
    nCen=rep(0,n)   
    Y=rep(0,n)      
    ll=rep(0,n)
    
    prod=1
    summa=0
    R0=n
    d=df$d
    t=df$t
    ties0=0
    nCen0=0
    ev0=0
    if(adjvar){
      Delta=df$Delta
      Q=df$Q
      Q2=df$Q2
    }
    
    for (i in 1:n){
      
      #conteggio ties
      if(i!=n){      
        if(t[i+1]==t[i]){ties0=ties0+1
        }else{           ties0=0}
        ties[i]=ties0
      }
      
      if(adjvar){
        Y0=(1-Delta[i])+Q[i]
        haz=(Y0*d[i])/R0
      }else{
        Y0=d[i]
        haz=Y0/R0
      }
      
      fatt=1-haz
      add=ifelse(haz>=(1-zn),NA,haz/(fatt*R0))
      
      adj=0
      if(adjvar) adj=ifelse(haz>=(1-zn),NA,((1-haz)**-2) * (1-1/R0) * ((Q[i]-Q2[i])/R0**2))
      
      prod=prod*fatt
      summa=summa+add+adj  
      
      Y[i]=Y0
      S[i]=prod
      R[i]=R0
      
      cv2[i]=ifelse(S[i]<=zn,NA,summa)
      #print(haz)
      ll[i]=ifelse(haz>=(1-zn) || haz<=zn
                   ,0,d[i]*log(haz)+(R0-d[i])*log(1-haz))
      
      DeltaR=1
      if(adjvar) DeltaR=1+(Y0-1)*d[i]
      R0=R0-DeltaR
      
      # ties and Cen counted
      if(i!=1){      
        if(t[i-1]==t[i]){ev0=ev0+Y0
        }else{           ev0=Y0}
        Y[i]=ev0
        if(d[i-1]==0){nCen0=nCen0+1
        }else{        nCen0=0}
        nCen[i]=nCen0
        R[i]=R[i]+ties[i-1] 
      }
    }
    cv=sqrt(cv2)
    sd_S=S*cv
    dfout=data.frame(df,R,Y,nCen,S,sd_S,cv,ties,ll,ll_sum=sum(ll))
    
    if(delcen) {
      cond=c(FALSE,dfout$ties[-n]==1 & dfout$d[-1]==0)
      dfout$d[cond]=1 
      dfout=dfout[dfout$d==1,]
    }    
    if(delties){dfout=dfout[dfout$ties==0,]} # elimina ties
    if(compfinpoint){  
      t.death.max=max(df[df$d==1,]$t) 
      t.lfu.res=df[df$t>t.death.max,]$t 
      if(delcen){
        if (length(t.lfu.res)>0){
          col=colnames(dfout)
          temp=data.frame(t(rep(NA,length(col))))
          colnames(temp)=col
          temp$t=max(t.lfu.res)
          temp$S=0
          dfout=rbind(dfout,temp)
        }
      }else{ 
        dfout[nrow(dfout),]$S=0 
      }
    }  
    
    dfout$p=-diff(c(1,dfout$S))  #compute probabilities
    
    row.names(dfout)=1:nrow(dfout)
    return(dfout)
  }
  
  ###############   comp_ciKM  function  #############
  # Compute confidence interval for KM estimators
  # input: df: data.frame(S,sd_S,cv)
  #        type.ci="log","linear", conf_lev=0.95
  # output: data.frame(ci_lower_S,ci_upper_S)
  comp_ciKM=function(df,
                     type.ci="log",
                     conf_lev= 0.95,
                     restrict=TRUE){
    if (conf_lev < 0 || conf_lev > 1) stop("confidence level must be between 0 and 1")
    
    zalpha <- qnorm(1-(1-conf_lev)/2)
    
    S=df$S
    cv=df$cv
    sd_S=df$sd_S
    R=df$R
    
    #log - greenwood
    if(type.ci=="log"){
      lower <- exp(log(S) - zalpha * cv)
      upper <- exp(log(S) + zalpha * cv)
    }
    if(type.ci=="linear"){
      lower <- S - zalpha * sd_S
      upper <- S + zalpha * sd_S
    }
    #peto
    if(type.ci=="peto"){
      peto.se <- cv/log(S)
      binom.se <- S * sqrt((1 - S)/R)
      lower <- S - zalpha * binom.se
      upper <- S + zalpha * binom.se
    }
    #loglog
    if(type.ci=="loglog"){
      peto.se <- cv/log(S)
      lower <- S**(exp(-zalpha * peto.se))
      upper <- S**(exp( zalpha * peto.se))
    }
    
    if (restrict){
      lower[lower < 0]=0
      upper[upper > 1]=1
    }
    out=data.frame(ci_lower_S=lower,
                   ci_upper_S=upper)
    return(out)
  }
  
  ###############   compute_q_ij  function #############
  # input: pdf: data.frame(t,p),
  #        theta
  # output: q_ij 
  compute_q_ij=function(pdf,theta){
    pdf=pdf[order(pdf$t),]
    nc=length(theta)
    nr=nrow(pdf)
    q_ij=matrix(0,nr,nc)
    q_compl=pdf$p
    for(j in 1:nc){
      cond=(pdf$t>theta[j])
      q_ij[cond,j]=q_compl[cond]
      q_ij[,j]=q_ij[,j]/sum(q_ij[,j])
    }
    return(q_ij)
  }
  
  #############   sample_cond_t function   ############
  # Draw a value of t from a conditional table given t>theta
  # input:  df: dataframe(t,p)
  #         theta: scalar
  # output: t_theta
  sample_cond_t=function(df,theta){
    df_res=df[df$t>theta,]
    if(nrow(df_res)>0){
      if(nrow(df_res)>1){t_theta=sample(df_res$t,size =1,replace = TRUE,prob = df_res$p)
      }else{             t_theta=df_res$t[1]}
    }else{
      t_theta=max(df$t)
    }
    return(t_theta)
  }
  
  ####  loess_h function    ########################
  # Local polynomial regression fitting  for 
  # incremental Nelson-Aalen cumulative hazard estimate  
  loess_h=function(df){
    outKM <- survival::survfit(Surv(df$t, df$d) ~ 1, ctype=1) 
    t <- outKM$time
    h <- diff(c(0,outKM$cumhaz))/diff(c(0,outKM$time))
    out=try(smooth <- stats::loess.smooth(t, h),silent=TRUE)
    if (!class(out)=="try-error"){
      
    }else{
      plus=0
      while(class(out)=="try-error" & plus<3){
        plus=plus+0.01
        out=try(smooth <- stats::loess.smooth(t, h+plus),silent=TRUE)
      }
    }
    hs=smooth$y
    ts=smooth$x
    return(data.frame(ts,hs))
  }
  
  ##############  alpha_loess  function    ############
  alpha_loess=function(df,theta){
    h1=loess_h(df)
    df0=df
    df0$d=1-df0$d
    h0=loess_h(df0)
    # ridefiniti se minori di zero
    if(sum(h1$hs<0)>0) h1[h1$hs<0,]$hs=0
    if(sum(h0$hs<0)>0) h0[h0$hs<0,]$hs=0
    # ridefiniti se maggiori di uno zero
    if(sum(h1$hs>1)>0) h1[h1$hs>1,]$hs=1
    if(sum(h0$hs>1)>0) h0[h0$hs>1,]$hs=1
    alpha=h1$hs/(h1$hs+h0$hs)
    # in caso di denominatore nullo ponpi a 0.5
    alpha[(h1$hs+h0$hs)==0]=0.5
    alpha_j=approx(h1$ts,alpha,theta)$y
    # ridefiniti se theta<min(h1$ts) e se theta>max(h1$ts)
    if(theta<min(h1$ts)) alpha_j=alpha[1]
    if(theta>max(h1$ts)) alpha_j=alpha[length(alpha)]
    return(alpha_j)
  }
  
  ##############    compute_alpha function   ##############
  compute_alpha=function(df_t,theta,type.alpha){
    nt=length(theta)
    if(type.alpha==4){
      alpha_j=sapply(1:nt,function(x) alpha_loess(df_t,theta[x]))
    }
    if(type.alpha==5){
      alpha_j=rep(0,nt)
    } 
    return(alpha_j)
  }
}
