#' Plots for betaclust object
#' @export
#' @param object betaclust object
#' @param what The different plots that can be obtained from the object (default="density") (what=c("density","uncertainty","InformationCriterion"))
#' @param plot_type The plot type to be displayed (default="ggplot")(plot_type="ggplot" or"plotly")
#' @param scale_param The axis that needs to be fixed or not for facet plot (default="free_y") (scales=c("free_y","free_x","free"))
#' @importFrom ggplot2 ggplot aes
#' @importFrom  plotly ggplotly

plot.betaclust <- function(object,what="density",
                                    plot_type="ggplot",scale_param="free_y")
{
  if(what=="density")
  {
    data<-object$best_model$data
    data_ggplot<-as.data.frame(data)
    data_ggplot$mem_final<-as.factor(data_ggplot$mem_final)
    colnames(data_ggplot)[length(data_ggplot)]<-"Cluster"
    if(object$optimal_model == "C.." || object$optimal_model == "CN.")
    {
      plot_density<-ggplot2::ggplot(data_ggplot,ggplot2::aes(x=data_ggplot[,1], fill=Cluster))+
        ggplot2::geom_density(alpha=0.6)+
        ggplot2::labs(x="Beta value", y="Density",
                      title=paste0("Density estimates for ",object$optimal_model,"  clustering solution"),
                      fill ="Cluster")
      #plot_density
      if(plot_type=="ggplot")
        plot_density
      else
        plotly::ggplotly(plot_density)


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
      color_length<-col_len-1
      colours<-scales::seq_gradient_pal(low="#FFC20A",high="#0C7BDC",space = "Lab")(1:color_length/color_length)

      plot_density<-ggplot2::ggplot(data_plot)+
        ggplot2::geom_density(aes(x=beta_value,color=Patient_Sample))+
        ggplot2::xlab("Beta Value")+
        ggplot2::ylab("Density")+
        ggplot2::scale_color_manual(values=colours)+
        ggplot2::facet_wrap(~Cluster,scales = scale_param
                            # ,labeller = ggplot2::labeller(Cluster= cluster_count_label)
        )+
        ggplot2::theme(axis.title.x = ggplot2::element_text(size=10),
                       axis.title.y = ggplot2::element_text(size=10)) +
        ggplot2::ggtitle("Density estimates for C.R clustering solution")

      #plot_density
      if(plot_type=="ggplot")
        plot_density
      else
        plotly::ggplotly(plot_density)


    }
  }else if(what=="uncertainty")
  {
    unc_df<-cbind(object$best_model$uncertainty,object$best_model$data[,"mem_final"])
    unc_df<-as.data.frame(unc_df)
    colnames(unc_df)<-c("Uncertainty","Cluster")
    unc_df$Cluster<-as.factor(unc_df$Cluster)
    unc_df_sorted<-unc_df[order(unc_df$Cluster),]
    max_unc=1-1/(length(object$best_model$tau))
    plot_uncertainty<-ggplot2::ggplot(unc_df_sorted, ggplot2::aes(x=Cluster, y=Uncertainty, color=Cluster)) +
      ggplot2::geom_boxplot()+ ggplot2::theme(legend.position = "none")+
      ggplot2::ggtitle("Boxplot for uncertainties in clustering solution")+
      ggplot2::coord_cartesian(ylim = c(0, 1))+
      ggplot2::geom_hline(yintercept=max_unc)+
      ggplot2::annotate(geom="text", label="maximum uncertainty",
               x=3, y=(max_unc+0.015), vjust=-1)

    if(plot_type=="ggplot")
      plot_uncertainty
    else
      plotly::ggplotly(plot_uncertainty)

  }else
  {
    if(length(object$ic_output)>1)
    {
      model_name<-as.vector(object$function_call$model_names[2:length(object$function_call$model_names)])
      model_name_wo_dot<-gsub("[[:punct:]]", "", model_name)
      ic_op<-object$ic_output
      ic_df<-do.call(rbind, Map(data.frame, A=ic_op, B=model_name_wo_dot))
      colnames(ic_df)<-c("IC_value","ModelName")
      ic_df2<-ic_df[order(ic_df$IC_value),]
      plot_ic<-ggplot2::ggplot(data=ic_df,ggplot2::aes(x=ModelName,y=IC_value,group=1))+
        ggplot2::geom_line()+
        ggplot2::ggtitle(paste0(wrapper_out2$information_criterion," Information Criterion Plot for optimal model selection"))+
        ggplot2::xlab("Model Name")+
        ggplot2::ylab("Information criterion value")

      if(plot_type=="ggplot")
        plot_ic
      else
        plotly::ggplotly(plot_ic)
    }else
      stop("Information criteria graph cannot be plotted for a single model.
      \n More than one model needs to be run to visualize this graph")

  }

}
