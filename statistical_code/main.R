run_protein_analysis <- function(input_csv_file, output_csv_directory) {
  # the R "glmnet" package needs to be installed,  https://cran.r-project.org/web/packages/glmnet/index.html
  library(glmnet); library(MASS)  
  
  # GNU Fortran needs also be installed, https://gcc.gnu.org/fortran/,  https://capsis.cirad.fr/capsis/documentation/mingw-installation
  # the file "nodelasso.dll" needs be stored first in disk, for example, the following show stored under Documents\protein folder
  dyn.load("C:\\Shabani\\Projects\\\tumor_latency\\statistical_code\\nodelasso.dll")
  
  nl2=100; lmax_node=10; lambds=lmax_node*(0.001)^(c(0:(nl2-1))/(nl2-1));  iters=100; ## creating one hundred numbers between 0 and 10. 
  
  # Read data file, the following example is for reading  data file “data_protein.csv”
  # In the data file, the first variable is the group indicator (e.g, 1: radiator; 0; Control), the rest variables are protein 
  datas=data.matrix(read.csv(file= input_csv_file, header = FALSE))
  
  cuts=0.00 
  
  y=datas[,1];      # Y is the response variable:  (e.g., 1: radiator; 0; Control)
  x1=datas[,-1];    # x1 is the matrix of proteins (e.g., if N=30, there are 1000 proteins, the matrix is 30*1000)
  nn=length(y); p1=length(x1[1,])
  
  
  fit1=glmnet(datas[,-1],datas[,1],family="binomial"); lglmnet=sort(fit1$lambda)
  lambdacv=cv.glmnet(datas[,-1],datas[,1],family="binomial",type.measure="class",lambda=lglmnet)
  lss=lambdacv$lambda.min
  
  beta_K=as.vector(coef(fit1,s=lss))  
  betaest=beta_K[-1]; index=which(abs(betaest)>cuts);  # Index is a vector show the positions of proteins, which have non-zero estimated regression coefficients 
  p2=length(index); 
  
  pii=exp(beta_K[1]+x1%*%beta_K[-1])/(1+exp(beta_K[1]+x1%*%beta_K[-1]))
  
  ## kernel gradient
  ker_grad=matrix(0,p1,p1)
  for (K in 1:nn){ker_grad=ker_grad+(x1[K,]*(y[K]-pii[K]))%*%(t(x1[K,]*(y[K]-pii[K])))}
  ker_grad=ker_grad/nn
  
  ## calculate gradient and hessisn
  grad=c(1:p1)*0; hessian=matrix(0,p1,p1)
  for (K in 1:p1){
    grad[K]=sum(x1[,K]*(y-pii))/nn
    for (L in K:p1){hessian[K,L]=(-1)*sum(x1[,K]*x1[,L]*pii*(1-pii))/nn}
  }
  grad=(-1)*grad; hessian=hessian+t(hessian); diag(hessian)=diag(hessian)/2;  hessian=(-1)*hessian
  
  chat=matrix(1,p1,p1); T2=chat*0;
  
  for (K in 1:p1){
    reco7=matrix(0,nl2,p1+1)
    for (L in 1:nl2){
      est=.Fortran("nodelasso",as.integer(p1),as.integer(K),as.double(hessian[K,-K]),hessian[-K,-K],as.double(c(1:(p1-1))*0),as.double(lambds[L]),as.integer(iters))
      gamma_K=est[[5]]
      reco7[L,1:(p1-1)]=gamma_K
      reco7[L,p1]=lambds[L]
      reco7[L,p1+1]=(-2)*(hessian[K,K]-2*(hessian[K,-K]%*%gamma_K)+t(gamma_K)%*%(hessian[-K,-K]%*%gamma_K))-length(gamma_K[which(gamma_K!=0)])*log(p1-1)
    }
    if(length(which(reco7[,p1+1]==max(reco7[,p1+1])))==1) {gamma_K=reco7[which(reco7[,p1+1]==max(reco7[,p1+1])),1:(p1-1)]
    }else{gamma_K=reco7[which(reco7[,p1+1]==max(reco7[,p1+1])),1:(p1-1)][1,]}
    chat[K,-K]=(-1)*gamma_K; T2[K,K]=1/(hessian[K,K]-hessian[K,-K]%*%gamma_K);
  }
  
  theta_l=T2%*%chat; biast=theta_l%*%grad; maxd=max(abs(theta_l%*%hessian-diag(1,p1))); dspn=betaest-biast;
  sigb=theta_l%*%((ker_grad)%*%t(theta_l)); sigma=diag(sigb)/nn
  
  estds=matrix(0,nrow=p2,ncol=7); estds[,1]=index;
  
  for (II in 1:p2){
    estds[II,2]=dspn[index[II]]; estds[II,3]=sqrt(sigma[index[II]]);  estds[II,4]=dspn[index[II]]-qnorm(0.975)*sqrt(sigma[index[II]]);  estds[II,5]=dspn[index[II]]+qnorm(0.975)*sqrt(sigma[index[II]]); estds[II,6]=(dspn[index[II]]/sqrt(sigma[index[II]]))^2
    estds[II,7]=2*pnorm(abs(dspn[index[II]]/sqrt(sigma[index[II]])),lower.tail = FALSE, log.p = FALSE)
  }  
  # Write the estds data frame to the output file
  write.csv(estds, file = output_csv_directory, row.names = FALSE)
}

# Usage example:

# List all CSV files in the directory where your input files are located
input_directory = "C:\\Shabani\\Projects\\tumor_latency\\data\\output_preprocess"
input_csv_files <- list.files(path = input_directory, pattern = ".csv", full.names = TRUE)
print(input_csv_files)


# Specify the directory where you want to save the output CSV files
output_csv_directory <- "C:\\Shabani\\Projects\\tumor_latency\\data\\output_r"

# Loop through each input CSV file and apply the function
for (input_file in input_csv_files) {
  ## print the links 
  print(input_file)
  
  ## Generate the output file name based on the input file name
  output_file <- file.path(output_csv_directory, paste0("output_", basename(input_file)))
  
  print(output_file)
  
  ## run the main function 
  run_protein_analysis(input_file, output_file)
  
  ## Print a message to indicate progress
  cat("Processed:", input_file, "\n")
}
