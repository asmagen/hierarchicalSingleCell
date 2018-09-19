#' return the pairwise Jaccard coefficients between two experiemnts
#' @param inputA a data frame with column names node (the name of the cluster) and
#' gene (all overexpressed gene of the node)
#' @param inputB same format as inputA
#' @return a list containing the Jaccard coefficients in long and wide form
get_jaccard<-function(inputA, inputB){
 splitA<-split(inputA,inputA$node)
 splitB<-split(inputB,inputB$node)
 num<-as.numeric(length(splitA)*length(splitB))
 jc<-data.frame(matrix(ncol=3,nrow=num))
 k=0
 for (i in seq(length(splitA))){
   v1<-splitA[[i]]$gene
   for (j in seq(length(splitB))){
     k=k+1
     print
     v2<-splitB[[j]]$gene
     jaccard_calc<-length(intersect(v1,v2))/length(union(v1,v2))
     jc[k,]<-c(names(splitA)[i],names(splitB)[j],jaccard_calc)
   }
 }
 jc$X3<-as.numeric(jc$X3)
 jc2<-matrix(jc$X3,nrow=length(unique(jc$X1)),ncol=length(unique(jc$X2)), byrow = TRUE)
 colnames(jc)<-c("node_exp1","node_exp2","j_indx")
 rownames(jc2)<-names(splitA)
 colnames(jc2)<-names(splitB)
 return(list(jc=jc,jc2=jc2))
}


