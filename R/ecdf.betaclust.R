globalVariables(c("beta_value","Patient_Sample"))
#' @title The empirical cumulative distribution function plot
#' @description An empirical cumulative distribution function (ECDF) plot for a \code{\link[betaclust:betaclust]{betaclust}} object.
#' @details This function plots the ECDF of the differentially methylated CpG sites identified using the K.R model for all patient samples.
#' The plot visualises the methylation state changes between the different DNA samples for each patient.
#' @export
#' @param x A dataframe containing methylation values of identified differentially methylated regions related to a gene.
#' Samples are grouped together in the dataframe such that the columns are ordered as Sample1_Patient1, Sample1_Patient2, Sample2_Patient1, Sample2_Patient2, etc.
#' @param R Number of tissue samples from which DNA methylation data are collected (default R = 2).
#' @param sample_name The order in which the samples are grouped in the dataframe. If no value is specified then default values of sample names, e.g. Sample 1, Sample 2, etc are used (default = NULL).
#' @param title The title that the user wants to display on the graph. The default is "NULL".
#' @return The ECDF plot for the selected CpG sites for all patients and their DNA samples.
#'
#' @seealso \code{\link{betaclust}}
#' @seealso \code{\link{beta_kr}}
#' @importFrom ggplot2 ggplot aes
#' @importFrom  plotly ggplotly
#' @importFrom scales seq_gradient_pal

ecdf.betaclust <- function(x, R=2, sample_name=NULL,title=NULL){

  if(is.null(title))
  {
    txt=""
  }else{
    txt=title
  }
  if(is.null(sample_name))
  {
    sample_name=sapply(1:R, function(x) paste0("Sample ",x))
  }
  col_len<-ncol(x)
  row_len<-nrow(x)
  col_names<-colnames(x)
  data_matrix<-as.matrix(x[,1:col_len])
  data_new<-as.vector(data_matrix)
  ecdf_df<-as.data.frame(data_new)
  Patient_sample=vector()
  Samples=vector()
  colnames(ecdf_df)[1]<-"beta_value"
  for(i in 1:(col_len))
  {
    temp=gsub("_"," ",col_names[i])
    ps_names<-rep(temp,row_len)
    Patient_sample<-c(Patient_sample,ps_names)
  }
  for(j in 1:R)
  {
    for(i in 1:(col_len/R))
    {
      samp_names<-rep(sample_name[j],row_len)
      Samples<-c(Samples,samp_names)
    }
  }

  data_plot<-cbind(ecdf_df,Patient_sample,Samples)
  data_plot<-as.data.frame(data_plot)
  colnames(data_plot)<-c("beta_value","Patient_Sample","Samples")

  data_plot$beta_value<-as.numeric(data_plot$beta_value)




   if(R<3)
   {
    patient_count<-col_len/R
    colours<-vector()
    color_length<-length(sample_name)
    colors<-scales::seq_gradient_pal(low="#FFC20A",high="#0C7BDC",space = "Lab")(1:color_length/color_length)
    colours<-sapply(1:color_length, function(x) rep(colors[x],times=patient_count))
    colours<-as.vector(colours)
    pecdf<-ggplot2::ggplot(data_plot, ggplot2::aes(beta_value, colour = Patient_Sample)) +
      ggplot2::stat_ecdf()+
      ggplot2::scale_color_manual("DNA Sample",values=colours)+
      ggplot2::ggtitle(txt)+
      #ggplot2::ggtitle("Empirical Cumulative Distribution Function")+
      ggplot2::xlab("Beta value")+
      ggplot2::ylab("F(Beta value)")
   }else if(R>2 & R<7){
     patient_count<-col_len/R
     colours<-vector()
     colours<-c(colours,rep("black",times=patient_count),rep("red",times=patient_count),rep("darkblue",times=patient_count),
                rep("darkgreen",times=patient_count),rep("yellow",times=patient_count),rep("brown",times=patient_count),
                rep("lavender",times=patient_count))
     colours<-colours[1:col_len]
     pecdf<-ggplot2::ggplot(data_plot, ggplot2::aes(beta_value, colour = Patient_Sample)) +
       ggplot2::stat_ecdf()+
       ggplot2::scale_color_manual(values=colours)+
       ggplot2::ggtitle(txt)+
       #ggplot2::ggtitle("Empirical Cumulative Distribution Function")+
       ggplot2::xlab("Beta value")+
       ggplot2::ylab("F(Beta value)")
  }
  else{
    pecdf<-ggplot2::ggplot(data_plot, ggplot2::aes(beta_value, colour = Patient_Sample)) +
      ggplot2::stat_ecdf()+
      ggplot2::labs(color="DNA Sample")+
      #ggplot2::scale_color_manual(values=colours)+
      #ggplot2::ggtitle("Empirical Cumulative Distribution Function")+
      ggplot2::ggtitle(txt)+
      ggplot2::xlab("Beta value")+
      ggplot2::ylab("F(Beta value)")
  }


  pecdf
}
