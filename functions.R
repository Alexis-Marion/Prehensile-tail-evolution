'%!in%' <- function(x,y)!('%in%'(x,y))

red<-rgb(238,99,99,alpha=0.8*255,maxColorValue=255)
blue<-rgb(100,149,237,alpha=0.8*255,maxColorValue=255)
green<-rgb(162,205,90,alpha=0.8*255,maxColorValue=255)


pGLS.plotGrade<-function (Yvar, Xvar, data, tree, model, group, ...) 
{
    dataTemp <- pruneSample(na.omit(data), tree, group)$data
    treeTemp <- pruneSample(dataTemp, tree, group)$tree
    Y <- dataTemp[, which(colnames(dataTemp) == paste(Yvar))]
    X <- dataTemp[, which(colnames(dataTemp) == paste(Xvar))]
    dataGLS <- as.data.frame(cbind(Y, X))
    rownames(dataGLS) <- rownames(dataTemp)
     switch(model, BM = {
        pGLSTemp <- phylolm(Y ~ X, data=dataGLS, phy=treeTemp,model = "BM")
        }, lambda = {
        pGLSTemp <- phylolm(Y ~ X, data=dataGLS, phy=treeTemp,model = "lambda")
        })
    a <- summary(pGLSTemp)$coefficients[1, 1]
    b <- summary(pGLSTemp)$coefficients[2, 1]
    lines(c(min(X), max(X)), c((a + b * min(X)), (a + b * max(X))), 
        ...)
    points(X, Y, ...)
}

getParameters<-function(data,tree,Group){
dataGroup<-pruneSample(data,tree,Group)$data
treeGroup<-pruneSample(data,tree,Group)$tree
      pglsModelGroup<-phylolm(Dependent~Independent,dataGroup,treeGroup,model="lambda")
      
        Slope<-summary(pglsModelGroup)$coefficients[2,1]; SlopeStderr<-summary(pglsModelGroup)$coefficients[2,2]; Intercept<-summary(pglsModelGroup)$coefficients[1,1]; InterceptStderr<-summary(pglsModelGroup)$coefficients[1,2]
        SlopeMIN<-Slope-qt(0.975,pglsModelGroup$n)*SlopeStderr
        SlopeMAX<-Slope+qt(0.975,pglsModelGroup$n)*SlopeStderr
        InterceptMIN<-Intercept-qt(0.975,pglsModelGroup$n)*InterceptStderr
        InterceptMAX<-Intercept+qt(0.975,pglsModelGroup$n)*InterceptStderr
    resultsSlopes<-rbind(Slope,SlopeMIN,SlopeMAX)
    resultsIntercepts<-rbind(Intercept,InterceptMIN,InterceptMAX)
    results<-rbind(resultsSlopes,resultsIntercepts)
      return(results)}

getGroupParameters<-function(data,tree,Group){
dataGroup<-pruneSample(data,tree,Group)$data
treeGroup<-pruneSample(data,tree,Group)$tree
      pglsModelGroup<-phylolm(Dependent~Independent,dataGroup,treeGroup,model="lambda")
      return(summary(pglsModelGroup)$coefficients)}

gls.Ymean<-function(Y,Sigma){
      n<-length(Y)
#Standardize sigma to reduce rounding errors
#      tr<-sum(diag(Sigma))
#      Sigma<-n*Sigma/tr
#Input
      invSigma<-solve(Sigma)
#pgls mean and variance
      X1<-rep(1,n)
      q<-2          # correct if multivariate!!
      C1<-solve(t(X1)%*%invSigma%*%X1)
      Y_PGLSmean<-C1%*%t(X1)%*%invSigma%*%Y
      Y_PGLSdeviations = Y - c(Y_PGLSmean)
      Y_PGLSvariance = (t(Y_PGLSdeviations)%*%invSigma%*%Y_PGLSdeviations)/(n-1)
      SE_Y_mean = sqrt(Y_PGLSvariance/n)
#Save model
      results<-cbind(Y_PGLSmean,SE_Y_mean,Y_PGLSvariance)
      colnames(results)<-c("Ymean","YSE","Y_PGLSvariance")
return(results)
}


getMeansDependent<-function(data,tree,Group){
      
data_Group<-pruneSample(data,tree,Group)$data
tree_Group<-pruneSample(data,tree,Group)$tree
Y<-data_Group$Dependent; Sigma<-vcv(tree_Group)
ResultsDependent<-gls.Ymean(Y,Sigma)
return(ResultsDependent)

}

getMeansIndependent<-function(data,tree,Group){
      
data_Group<-pruneSample(data,tree,Group)$data
tree_Group<-pruneSample(data,tree,Group)$tree
Y<-data_Group$Independent; Sigma<-vcv(tree_Group)
ResultsIndependent<-gls.Ymean(Y,Sigma)
return(ResultsIndependent)
}


getpGLSgroup<-function(targetGroup){
  dataTemp<-pruneSample(data,tree,targetGroup)$data
  treeTemp<-pruneSample(data,tree,targetGroup)$tree
  pGLSTemp<-summary(phylolm(Dependent~Independent,dataTemp,treeTemp,model="lambda"))
		        Slope<-pGLSTemp$coefficients[2,1]; SlopeStderr<-pGLSTemp$coefficients[2,2]; Intercept<-pGLSTemp$coefficients[1,1]; InterceptStderr<-pGLSTemp$coefficients[1,2]
		        N<-length(treeTemp$tip.label)-2
		        SlopeMIN<-Slope-qt(0.975,N)*SlopeStderr
		        SlopeMAX<-Slope+qt(0.975,N)*SlopeStderr
		        InterceptMIN<-Intercept-qt(0.975,N)*InterceptStderr
		        InterceptMAX<-Intercept+qt(0.975,N)*InterceptStderr
		    resultsSlopes<-cbind(Slope,SlopeMIN,SlopeMAX)
		    resultsIntercepts<-cbind(Intercept,InterceptMIN,InterceptMAX)
		    results<-cbind(resultsSlopes,resultsIntercepts)
  return(results)
}
