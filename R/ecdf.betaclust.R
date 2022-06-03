#' @title The empirical cumulative distribution function
#' @description Empirical Cumulative Distribution Function plot for betaclust object
#' @export
#' @param x Methylation values of Identified Differentially methylated regions related to a gene.
#' Group each sample together in the dataframe such that the columns are ordered as --> Sample1_P1, Sample1_P2, Sample2_P1, Sample2_P2
#' @param samples number of tissue samples from where DNA methylation data is collected (default samples=2)
#' @param sample_name The order in which the samples are grouped in the dataframe (default = c("Sample 1","Sample 2"))
#' @return The ecdf plot for the selected CpG sites for all patients and samples.
#'
#' @importFrom ggplot2 ggplot aes
#' @importFrom  plotly ggplotly

ecdf.betaclust <- function(x, samples=2, sample_name = c("Sample 1","Sample 2")){

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
    ps_names<-rep(col_names[i],row_len)
    Patient_sample<-c(Patient_sample,ps_names)
  }
  for(j in 1:samples)
  {
    for(i in 1:(col_len/samples))
    {
      samp_names<-rep(sample_name[j],row_len)
      Samples<-c(Samples,samp_names)
    }
  }

  data_plot<-cbind(ecdf_df,Patient_sample,Samples)
  data_plot<-as.data.frame(data_plot)
  colnames(data_plot)<-c("beta_value","Patient_Sample","Samples")

  data_plot$beta_value<-as.numeric(data_plot$beta_value)




  if(samples<7)
  {
    patient_count<-col_len/samples
    colours<-vector()
    colours<-c(colours,rep("black",times=patient_count),rep("red",times=patient_count),rep("darkblue",times=patient_count),
               rep("darkgreen",times=patient_count),rep("yellow",times=patient_count),rep("brown",times=patient_count),
               rep("lavender",times=patient_count))
    colours<-colours[1:col_len]
    pecdf<-ggplot2::ggplot(data_plot, ggplot2::aes(beta_value, colour = Patient_Sample)) +
      ggplot2::stat_ecdf()+
      ggplot2::scale_color_manual(values=colours)+
      ggplot2::ggtitle("Empirical Cumulative Distribution Function")+
      ggplot2::xlab("Beta value")+
      ggplot2::ylab("F(Beta value)")
  }else{
    pecdf<-ggplot2::ggplot(data_plot, ggplot2::aes(beta_value, colour = Patient_Sample)) +
      ggplot2::stat_ecdf()+
      #ggplot2::scale_color_manual(values=colours)+
      ggplot2::ggtitle("Empirical Cumulative Distribution Function")+
      ggplot2::xlab("Beta value")+
      ggplot2::ylab("F(Beta value)")
  }


  pecdf
}
