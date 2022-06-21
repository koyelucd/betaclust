#' @title Plots for visualizing the betaclust class object
#' @description This function helps visualise the clustering solution by plotting the density estimates, the uncertainty and the information criterion.
#' @details The density estimates under the optimal clustering solution by specifying what = "density" in the function. Interactive plots can also be produced using
#'  plot_type = "plotly". The uncertainty in the clustering solution can be plotted using what="uncertainty".
#' The information criterion values for all fitted models can be plotted using what = "InformationCriterion".
#' @export
#' @seealso \code{\link{betaclust}}
#' @param object A betaclust object.
#' @param what The different plots that can be obtained are either "density","uncertainty" or "InformationCriterion". (default="density").
#' @param plot_type The plot type to be displayed are either "ggplot" or "plotly". (default="ggplot").
#' @param title The title that the user wants to display on the graph. If no title is to be displayed the default is "NULL" value.
#' @param scale_param The axis that needs to be fixed for density estimates plot for visualizing the C.R clustering solution are either "free_y","free_x" or "free". (default = "free_y").
#' @importFrom ggplot2 ggplot aes
#' @importFrom  plotly ggplotly

plot.betaclust <- function(object,what="density",
                           plot_type="ggplot",title=NULL,scale_param="free_y")
{

  if(is.null(title))
  {
    txt=""
  }else{
    txt=title
  }
  if(what == "density")
  {
    data<-object$optimal_model_results$data
    data_ggplot<-as.data.frame(data)
    data_ggplot$mem_final<-as.factor(data_ggplot$mem_final)
    colnames(data_ggplot)[length(data_ggplot)]<-"Cluster"

    if(object$optimal_model == "C.." || object$optimal_model == "CN.")
    {
      plot_graph<-ggplot2::ggplot(data_ggplot,ggplot2::aes(x=data_ggplot[,1], fill=Cluster))+
        ggplot2::geom_density(alpha=0.6)+
        ggplot2::labs(x="Beta value", y="Density",
                      #title=paste0("Density estimates for ",object$optimal_model,"  clustering solution"),
                      title=txt,
                      fill ="Cluster")

      p.data <- ggplot2::ggplot_build(plot_graph)$data[[1]]


      p.text <- lapply(split(p.data, f = p.data$group), function(df){
        df[which.max(df$scaled), ]
      })
      p.text <- do.call(rbind, p.text)
      p.text$prop=p.text$n/(sum(p.text$n))


      plot_graph<-plot_graph + ggplot2::annotate('text', x = p.text$x,
                                                 y = p.text$y,
                                                 label = sprintf('%.3f',
                                                                 p.text$prop),
                                                 vjust = 0)


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

      plot_graph<-ggplot2::ggplot(data_plot)+
        ggplot2::geom_density(aes(x=beta_value,color=Patient_Sample))+
        ggplot2::xlab("Beta Value")+
        ggplot2::ylab("Density")+
        ggplot2::scale_color_manual(values=colours)+
        ggplot2::facet_wrap(~Cluster,scales = scale_param
                            # ,labeller = ggplot2::labeller(Cluster= cluster_count_label)
        )+
        ggplot2::theme(axis.title.x = ggplot2::element_text(size=10),
                       axis.title.y = ggplot2::element_text(size=10)) +
        ggplot2::ggtitle(txt)
        #ggplot2::ggtitle("Density estimates for C.R clustering solution")

      f_labels<-data.frame(Cluster=seq(1,length(object$optimal_model_results$cluster_size),by=1),label=as.vector(round((object$optimal_model_results$cluster_size/object$C),3)))
      plot_graph<-plot_graph+
        ggplot2::geom_text(x = 0.2, y = 1, ggplot2::aes(label = label), data = f_labels)

    }


  }
  if(what == "uncertainty")
  {
    labels<-c(max_uc="Maximum uncertainty")
    unc_df<-cbind(object$optimal_model_results$uncertainty,object$optimal_model_results$data[,"mem_final"])
    unc_df<-as.data.frame(unc_df)
    colnames(unc_df)<-c("Uncertainty","Cluster")
    unc_df$Cluster<-as.factor(unc_df$Cluster)
    unc_df_sorted<-unc_df[order(unc_df$Cluster),]
    max_unc=1-1/(length(object$optimal_model_results$tau))
    h=max_unc+0.015
    max_uncertainty <- data.frame(yintercept=h, max_uncertainty=factor(h))
    plot_graph<-ggplot2::ggplot(unc_df_sorted, ggplot2::aes(x=Cluster, y=Uncertainty, color=Cluster)) +
      ggplot2::geom_boxplot()+
      ggplot2::theme(axis.text.x=ggplot2::element_blank(),
                     axis.ticks.x=ggplot2::element_blank()
      )+
      #ggplot2::theme(legend.position = "none")+
      ggplot2::ggtitle(txt)+
      #ggplot2::ggtitle("Boxplot for uncertainties in clustering solution")+
      ggplot2::coord_cartesian(ylim = c(0, 1))+
      ggplot2::geom_hline(ggplot2::aes(yintercept=max_unc, linetype="Maximum uncertainty"),color="black")+
      ggplot2::scale_linetype_manual(name="", values = 2,guide = ggplot2::guide_legend(override.aes = list(color = "black")))



  }
  if(what == "InformationCriterion")
  {
    if(length(object$ic_output)>1)
    {
      model_name<-as.vector(object$function_call$model_names[2:length(object$function_call$model_names)])
      model_name_wo_dot<-gsub("[[:punct:]]", "", model_name)
      for(i in 1:length(model_name_wo_dot)){
        if(model_name_wo_dot[i]=="C")
        {
          model_name_wo_dot[i]="C.."
        }else if(model_name_wo_dot[i]=="CN")
        {
          model_name_wo_dot[i]="CN."
        }else if(model_name_wo_dot[i]=="CR")
        {
          model_name_wo_dot[i]="C.R"
        }
      }
      ic_op<-object$ic_output
      ic_df<-do.call(rbind, Map(data.frame, A=ic_op, B=model_name_wo_dot))
      colnames(ic_df)<-c("IC_value","ModelName")
      ic_df2<-ic_df[order(-ic_df$IC_value),]
      ic_df2$ModelName<-factor(ic_df2$ModelName,levels = ic_df2$ModelName)
      plot_graph<-ggplot2::ggplot(data=ic_df2,ggplot2::aes(x=ModelName,y=IC_value))+
        #ggplot2::geom_line()+
        ggplot2::geom_bar(stat="identity",fill="lightgreen")+
        ggplot2::geom_text(ggplot2::aes(label = sprintf("%.1f", IC_value), y= IC_value),  vjust = -1)+
        ggplot2::guides(fill="none")+
        ggplot2::ggtitle(txt)+
        #ggplot2::ggtitle(paste0(object$information_criterion," Information Criterion Plot for optimal model selection"))+
        ggplot2::xlab("Model Name")+
        ggplot2::ylab("Information criterion value")


    }else
      stop("Information criteria graph cannot be plotted for a single model.
          \n More than one model needs to be run to visualize this graph")

  }


  if(plot_type=="ggplot")
    plot_graph
  else
    plotly::ggplotly(plot_graph)

}
