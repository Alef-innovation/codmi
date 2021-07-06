if(!require(survival)) stop('Please run install.packages("survival")')
if(!require(stats)) stop('Please run install.packages("stats")')

source("codmi_functions.R")

##############################################################
##########          input  data                   ############

ncog = read.table("https://web.stanford.edu/~hastie/CASI_files/DATA/ncog.txt",sep=" ",header=TRUE)

input_data="ncogA"
if(input_data=="ncogA"){
  data.basis<- subset(ncog, arm == c('A'), select = c(t, d))
  t.covid <- c(250, 500, 750, 1000, 1250)
}
if(input_data=="ncogB"){
  data.basis<- subset(ncog, arm == c('B'), select = c(t, d))
  t.covid <- c(400, 800, 1200, 1600, 2000)
}
if(input_data=="Toy"){
  data.basis=data.frame(
    t = c(7, 34, 42, 63, 63, 74, 74, 90, 93, 100), 
    d = c(1,  1,  0,  1,  1,  1,  0,  1,  0,   0))
  t.covid <- c(50,91)
}

# Kaplan Meier
pdf.basis=KM_AV(data.basis)

#1) codmi 
out1 = codmi(data.basis,t.covid,
             epsilon=0.1)
out1$converge
out1$mean_iter
out1$table_fin
out1$theta_cen
out1$mean_fin
out1$pdf_fin

#2) codmi with cenadj=TRUE
out2 = codmi(data.basis,t.covid,
             epsilon=0.1,
             cenadj = TRUE)
out2$converge
out2$mean_iter
out2$table_fin
out2$theta_cen
out2$mean_fin
out2$pdf_fin

#3) codmi with exogenous alpha
out3 = codmi(data.basis,t.covid,
             cenadj = TRUE,epsilon=0.1,
             exo_alpha = c(0,0,1,0.4,0.7))
out3$converge
out3$mean_iter
out3$table_fin
out3$theta_cen
out3$mean_fin
out3$pdf_fin

