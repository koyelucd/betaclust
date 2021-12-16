#' Plots for betaclust object
#' @export
#' @param object betaclust object
#' @param what The different plots that can be obtained from the object (default="density") (what=c("density","uncertainty","InformationCriterion"))
#' @param plot_type The plot type to be displayed (default="ggplot")(plot_type=c("ggplot","plotly"))
#' @importFrom ggplot2 ggplot aes

plot.betaclust <- function(object,what="density",
                                    plot_type="ggplot")
{

  data<-object$best_model$data
  data_ggplot<-as.data.frame(data)
  data_ggplot$mem_final<-as.factor(data_ggplot$mem_final)
  colnames(data_ggplot)[length(data_ggplot)]<-"Cluster"
  if(object$optimal_model == "C.." || object$optimal_model == "CN.")
  {
    plot_density<-ggplot2::ggplot(data_ggplot,ggplot2::aes(x=P_Sam1_1, fill=Cluster))+
            ggplot2::geom_density(alpha=0.6)+
      ggplot2::labs(x="Beta value", y="Density",
                    title=paste0("Density estimates for ",object$optimal_model,"  clustering solution"),
                    fill ="Cluster")
    plot_density
  }else
  {
    cols=ncol(data_ggplot)
    rows=nrow(data_ggplot)
    data_matrix<-as.matrix(data_ggplot[,1:(cols-1)])
    data_new<-as.vector(data_matrix)
    col_names<-colnames(data_ggplot)
    col_len<-length(col_names)
    Cluster<-vector()
    Patient_sample=vector()
    for(i in 1:(col_len-1))
    {
      ps_names<-rep(col_names[i],rows)
      Patient_sample<-c(Patient_sample,ps_names)
      Cluster<-c(Cluster,data_ggplot[,cols])
    }
    data_plot<-cbind(data_new,Cluster,Patient_sample)
    data_plot<-as.data.frame(data_plot)
    colnames(data_plot)<-c("beta_value","Cluster","Patient_Sample")
    data_plot$beta_value<-as.numeric(data_plot$beta_value)
    # cluster_count<-length(table(data_ggplot[,cols]))
    # cluster_count_label<-vector()
    # for(j in 1:cluster_count)
    # {
    #   text_label=paste0(j,"=", "Cluster ",j)
    #   cluster_count_label<-c(cluster_count_label,text_label)
    # }
    color_length<-col_len-1
    colours<-scales::seq_gradient_pal(low="#FFC20A",high="#0C7BDC",space = "Lab")(1:color_length/color_length)

    plot_density<-ggplot2::ggplot(data_plot)+
      ggplot2::geom_density(aes(x=beta_value,color=Patient_Sample))+
      ggplot2::xlab("Beta Value")+
      ggplot2::ylab("Density")+
      ggplot2::scale_color_manual(values=colours)+
      ggplot2::facet_wrap(~Cluster ,scales = "free"
                          # ,labeller = ggplot2::labeller(Cluster= cluster_count_label)
                          )+
      ggplot2::theme(axis.title.x = ggplot2::element_text(size=10),
            axis.title.y = ggplot2::element_text(size=10)) +
      ggplot2::ggtitle("Density estimates for C.R clustering solution")

    plot_density


  }


 # invisible()
}
