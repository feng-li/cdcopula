library("copula")
#try to combine them with varied copula
t.cop = tCopula(param=0.5, dim=2, df=4)
norm.cop = normalCopula(param=0.5, dim=2)
clayton.cop = archmCopula("clayton", param =0.5, dim = 2)
frank.cop = archmCopula("frank", param =0.5, dim = 2)
amh.cop = archmCopula("amh", param =0.5, dim = 2)
gumbel.cop = archmCopula("gumbel", param =5, dim = 2)
joe.cop = archmCopula("joe", param =5, dim = 2)

#reserve room
X = Y = matrix(0,1000,100)
VARX = VARY = matrix(0,1,100)

VAR = matrix(0,1,100)
VAR2 = matrix(0,1,100)
VAR3 = matrix(0,1,100)
VAR4 = matrix(0,1,100)

toplot = matrix(0,30,100)

#gengerate 100 groups and calculate their VaR
for(i in 1:100)
{
  set.seed(i)                   
  X[1:1000,i]=rt(1000,2)
  VARX[,i] = qt(0.95,2)     #calculate VaR of X to get the directional variable
  set.seed(i+1)
  Y[1:1000,i]=rt(1000,10)   
  VARY[,i] = qt(0.95,10)    #calculate VaR of Y to get the directional variable
}

#calculate the VaR of profolio
for(j in 1:100)
{
  mydata = matrix(0,1000,6)
  mydata[,1]=X[,j]
  mydata[,2]=Y[,j]
  dX=dt(X[,j],df=4)
  dY=dt(Y[,j],df=4)
  
  #use different copula
  for(i in 1:1000)
  {
    mydata[i,3] = pCopula(c(dX[i],dY[i]), t.cop)
  }
  
  for(i in 1:1000)
  {
    mydata[i,4] = pCopula(c(dX[i],dY[i]), norm.cop)
  }
  
  for(i in 1:1000)
  {
    mydata[i,5] = pCopula(c(dX[i],dY[i]), clayton.cop)
  }
  
  
  for(i in 1:1000)
  {
    mydata[i,6] = pCopula(c(dX[i],dY[i]), gumbel.cop)
  }
  
  #dangerous points
  outdata4 = mydata[which(mydata[,6]<0.01),]
  outdata = mydata[which(mydata[,3]<0.01),]
  outdata3 = mydata[which(mydata[,5]<0.01),]
  outdata2 = mydata[which(mydata[,4]<0.01),]

  #only hold the points with VaR<0
  okdata = outdata[which((VARX[,j]/(VARX[,j]+VARY[,j]))*outdata[,1]+(VARY[,j]/(VARX[,j]+VARY[,j]))*outdata[,2]<0),]
  
  #get the minimum distance of these points and the directional line
  abs=min(abs(okdata[,2]-(VARY[,j]/VARX[,j])*okdata[,1]))
  c=okdata[which(abs(okdata[,2]-(VARY[,j]/VARX[,j])*okdata[,1])==abs),]
  VAR[,j]=(VARX[,j]/(VARX[,j]+VARY[,j]))*c[1]+(VARY[,j]/(VARX[,j]+VARY[,j]))*c[2]
  VAR[,j]
  
  okdata2 = outdata2[which((VARX[,j]/(VARX[,j]+VARY[,j]))*outdata2[,1]+(VARY[,j]/(VARX[,j]+VARY[,j]))*outdata2[,2]<0),]
  
  abs2=min(abs(okdata2[,2]-(VARY[,j]/VARX[,j])*okdata2[,1]))
  c2=okdata2[which(abs(okdata2[,2]-(VARY[,j]/VARX[,j])*okdata2[,1])==abs2),]
  VAR2[,j]=(VARX[,j]/(VARX[,j]+VARY[,j]))*c2[1]+(VARY[,j]/(VARX[,j]+VARY[,j]))*c2[2]
  VAR2[,j]
  
  okdata3 = outdata3[which((VARX[,j]/(VARX[,j]+VARY[,j]))*outdata3[,1]+(VARY[,j]/(VARX[,j]+VARY[,j]))*outdata3[,2]<0),]
  
  abs3=min(abs(okdata3[,2]-(VARY[,j]/VARX[,j])*okdata3[,1]))
  c3=okdata3[which(abs(okdata3[,2]-(VARY[,j]/VARX[,j])*okdata3[,1])==abs3),]
  VAR3[,j]=(VARX[,j]/(VARX[,j]+VARY[,j]))*c3[1]+(VARY[,j]/(VARX[,j]+VARY[,j]))*c3[2]
  VAR3[,j]
  
  okdata4 = outdata4[which((VARX[,j]/(VARX[,j]+VARY[,j]))*outdata4[,1]+(VARY[,j]/(VARX[,j]+VARY[,j]))*outdata4[,2]<0),]
  
  abs4=min(abs(okdata4[,2]-(VARY[,j]/VARX[,j])*okdata4[,1]))
  c4=okdata4[which(abs(okdata4[,2]-(VARY[,j]/VARX[,j])*okdata4[,1])==abs4),]
  VAR4[,j]=(VARX[,j]/(VARX[,j]+VARY[,j]))*c4[1]+(VARY[,j]/(VARX[,j]+VARY[,j]))*c4[2]
  VAR4[,j]
  
  #take 30 points to plot
  toplot[,j] =sample( ((VARX[,j]/(VARX[,j]+VARY[,j]))*mydata[,1]+(VARY[,j]/(VARX[,j]+VARY[,j]))*mydata[,2]), 30)
  
}

plot(VAR[,1:100],type="l",ylim=c(-10,10))
lines(VAR2[,1:100],type="l",ylim=c(-10,10),col="red")
lines(VAR3[,1:100],type="l",ylim=c(-10,10),col="gold")
lines(VAR4[,1:100],type="l",ylim=c(-10,10),col="blue")


for(i in 1:100)
{
  a=rep(i,30)
  points(a,toplot[,i],col="green")
}






