globalVariables(c("Patient_Sample","label","Uncertainty","ModelName","IC_value"))
#' @title Plots for visualizing the betaclust class object
#' @description Visualise a \code{\link[betaclust:betaclust]{betaclust}} clustering solution by plotting the fitted and kernel density estimates, the uncertainty and the information criterion.
#' @details The fitted density estimates can be visualized under the optimal clustering solution by specifying what = "fitted density" and kernel density estimates under the optimal clustering solution by specifying what = "kernel density".
#'
#' The threshold inferred can be visualized by specifying threshold = TRUE.
#' The KN. model calculates different pairs of threshold points for each patient as the shape parameters are allowed to vary for each patient.
#' So the patient for whom the threshold needs to be displayed can be specified by inputting the column number representing the patient in the patient-wise ordered dataset in the parameter patient_number.
#'
#' Interactive plots can also be produced using plot_type = "plotly". The uncertainty in the clustering solution can be plotted using what = "uncertainty".
#' The information criterion values for all fitted models can be plotted using what = "information criterion".
#' @rdname plot.betaclust
#' @export
#' @examples
#' \donttest{
#' my.seed <- 190
#' M <- 3
#' N <- 4
#' R <- 2
#' data_output <- betaclust(pca.methylation.data[1:100,2:9], M, N, R,
#'             model_names = c("K..","KN.","K.R"), model_selection = "BIC",
#'             parallel_process = FALSE, seed = my.seed)
#' plot(data_output, what = "fitted density", plot_type = "ggplot",
#'      data=pca.methylation.data[1:100,2:9])
#' plot(data_output, what = "uncertainty", plot_type = "ggplot")
#' plot(data_output, what = "information criterion", plot_type = "ggplot")
#' }
#' @seealso \code{\link{betaclust}}
#' @param x A \code{\link[betaclust:betaclust]{betaclust}} object.
#' @param what The different plots that can be obtained are either "fitted density","kernel density","uncertainty" or "information criterion" (default = "fitted density").
#' @param plot_type The plot type to be displayed are either "ggplot" or "plotly" (default = "ggplot").
#' @param data A dataframe of dimension \eqn{C \times NR} containing methylation values for \eqn{C} CpG sites from \eqn{R} samples collected from \eqn{N} patients which was passed as an argument to the \code{\link[betaclust:betaclust]{betaclust}} function.
#' The data is not required as an input when generating "uncertainty" or "information criterion" plots and the default has been set as "NULL". The data needs to be passed as an argument to this function when generating either "fitted density" or "kernel density" plots.
#' @param sample_name The names of DNA samples in the dataset analysed by the K.R model. If no value is passed then default values of sample names, e.g. Sample 1, Sample 2, etc are used as legend text (default = NULL).
#' @param title The title that the user wants to display. If no title is to be displayed the default is "NULL".
#' @param patient_number The column number representing the patient in the patient-wise ordered dataset selected for visualizing the clustering solution of the K.. or KN. model (default = 1).
#' @param threshold The "TRUE" option displays the threshold points in the graph for the K.. and the KN. model (default = "FALSE").
#' @param scale_param The position scales can be fixed or allowed to vary between different panels generated for the density estimate plots for visualizing the K.R clustering solution. Options are "fixed", "free_y","free_x" or "free" (default = "free_y"). The option "fixed" results in the x and y scales being fixed across all panels, "free" varies the x and y scales across the panels, "free_x" fixes the y scale and lets the x scale vary across all panels and "free_y" fixes the x scale and lets the y scale vary across all panels.
#' @param ... Other graphics parameters.
#'
#' @return This function displays the following plots as requested by the user:
#' \itemize{
#' \item fitted density estimates - Plot showing the fitted density estimates of the clustering solution under the optimal model selected.
#' \item kernel density estimates - Plot showing the kernel density estimates of the clustering solution under the optimal model selected.
#' \item uncertainty -  A boxplot showing the uncertainties in the optimal clustering solution.
#' \item information criterion - Plot showing the information criterion values for all models fitted to support the selection of the optimal model.
#' }
#'
#'
#' @importFrom ggplot2 ggplot aes
#' @importFrom  plotly ggplotly
#' @importFrom stats C
#' @importFrom scales seq_gradient_pal

plot.betaclust <- function(x,what = "fitted density",
                           plot_type = "ggplot",data = NULL, sample_name = NULL,
                           title = NULL,patient_number = 1,
                           threshold = FALSE,scale_param = "free_y",...)
{
  object <- x
  pn=patient_number

  if(is.null(title))
  {
    txt=""
  }else{
    txt=title
  }
  if(is.null(sample_name))
  {
    R=object$R
    sample_name=sapply(1:R, function(x) paste0("Sample ",x))
  }
  if(what == "kernel density")
  {
    if(is.null(data)){
      plot_graph=NULL
      warning("data argument cannot be NULL for generation of kernel density plot.", call. = FALSE)
    }else
    {


      if(object$optimal_model == "K.." || object$optimal_model == "KN.")
      {
        #data<-object$optimal_model_results$data
        column_len=ncol(data)
        if(column_len>object$N)
        {
          call_data=data[,1:object$N]
        }else{
          call_data=data
        }
        data_ggplot<-as.data.frame(call_data)
        #data_ggplot$mem_final<-as.factor(data_ggplot$mem_final)
        data_ggplot$mem_final<-as.factor(object$optimal_model_results$classification)
        colnames(data_ggplot)[length(data_ggplot)]<-"Cluster"
        plot_graph<-ggplot2::ggplot(data_ggplot,ggplot2::aes(x=data_ggplot[,pn], fill=Cluster))+
          ggplot2::geom_density(alpha=0.6)+
          ggplot2::labs(x="Beta value", y="Density",
                        #title=paste0("Density estimates for ",object$optimal_model,"  clustering solution"),
                        title=txt,
                        fill ="Cluster")
        if(object$K==3)
        {
          colours<-c("chartreuse3","magenta","cyan3")
          plot_graph<-plot_graph+
            ggplot2::scale_fill_manual(values=colours)
        }

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
                                                   colour=p.text$fill,
                                                   vjust = 0)


      }else
      {
        column_len=ncol(data)
        if(column_len==(object$N*object$R))
        {
          call_data=data
        }else if(column_len>(object$N*object$R)){
          call_data=data[,1:(object$N*object$R)]
        }else{
          call_data=NULL
        }
        #data<-object$optimal_model_results$data
        if(is.null(call_data))
        {
          plot_graph=NULL
          warning("K.R considers only multiple samples. Please pass DNA methylation values for more than 1 sample.", call. = FALSE)
        }
        else{
        data_ggplot<-as.data.frame(call_data)
        #data_ggplot$mem_final<-as.factor(data_ggplot$mem_final)
        data_ggplot$mem_final<-as.factor(object$optimal_model_results$classification)
        colnames(data_ggplot)[length(data_ggplot)]<-"Cluster"
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
          temp=gsub("_"," ",col_names[i])
          ps_names<-rep(temp,rows)
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
          ggplot2::scale_color_manual("DNA sample",values=colours)+
          ggplot2::facet_wrap(~Cluster,scales = scale_param
                              # ,labeller = ggplot2::labeller(Cluster= cluster_count_label)
          )+
          ggplot2::theme(axis.title.x = ggplot2::element_text(size=10),
                         axis.title.y = ggplot2::element_text(size=10)) +
          ggplot2::ggtitle(txt)
        #ggplot2::ggtitle("Density estimates for K.R clustering solution")

        f_labels<-data.frame(Cluster=seq(1,length(object$optimal_model_results$cluster_size),by=1),label=as.vector(round((object$optimal_model_results$cluster_size/object$C),3)))
        plot_graph<-plot_graph+
          ggplot2::geom_text(x = 0.2, y = 1, ggplot2::aes(label = label), data = f_labels)

        }
      }

    }


  }
  if(what=="fitted density")
  {
    if(is.null(data)){
      plot_graph=NULL
      warning("data argument cannot be NULL for generation of fitted density plot.", call. = FALSE)
    }else
    {
      if(object$optimal_model == "K.." || object$optimal_model == "KN.")
      {

        #data_x=sort(object$optimal_model_results$data[,pn])
        data_x=sort(data[,pn])
        K=object$K
        prop=object$optimal_model_results$tau
        data_th_plot<-matrix(data=NA,nrow=1,ncol=3)
        data_th_plot<-as.data.frame(data_th_plot)
        colnames(data_th_plot)<-c("beta_value","density","Cluster")
        alpha=object$optimal_model_results$alpha
        delta=object$optimal_model_results$delta
        if(object$optimal_model == "K..")
        {
          for(i in 1:K)
          {
            beta_value=data_x
            density=prop[i]*stats::dbeta(data_x,alpha[i],delta[i])
            Cluster<-rep(i,length(data_x))
            temp<-cbind(beta_value,density,Cluster)
            data_th_plot<-rbind(data_th_plot,temp)
          }
        }else if(object$optimal_model == "KN.")
        {
          for(i in 1:K)
          {
            beta_value=data_x
            density=prop[i]*stats::dbeta(data_x,alpha[i,pn],delta[i,pn])
            Cluster<-rep(i,length(data_x))
            temp<-cbind(beta_value,density,Cluster)
            data_th_plot<-rbind(data_th_plot,temp)
          }
        }
        data_th_plot<-as.data.frame(data_th_plot)
        data_th_plot<-data_th_plot[-1,]
        data_th_plot$Cluster<-as.factor(data_th_plot$Cluster)
        #txt=""
        #colours<-c("chartreuse3","magenta","cyan3")
        plot_graph<-ggplot2::ggplot(data_th_plot)+
          ggplot2::geom_line(ggplot2::aes(beta_value,density,color=Cluster),
                             linetype = "solid")+
          ggplot2::labs(x="Beta value",y="Density",title=txt,
                        color ="Cluster")
        # +
        #   ggplot2::scale_color_manual(values=colours)
        if(K==3)
        {
          colours<-c("chartreuse3","magenta","cyan3")
          plot_graph<-plot_graph+
            ggplot2::scale_color_manual(values=colours)
        }
        p.data <- ggplot2::ggplot_build(plot_graph)$data[[1]]


        p.text <- lapply(split(p.data, f = p.data$group), function(df){
          df[which.max(df$y), ]
        })
        p.text <- do.call(rbind, p.text)
        p.text$prop=prop
        plot_graph<-plot_graph + ggplot2::annotate('text', x = p.text$x,
                                                   #y = p.text$y,
                                                   y=0.2,
                                                   label = sprintf('%.3f',
                                                                   p.text$prop),
                                                   colour=p.text$colour,
                                                   vjust = 0)
      }
      else{
        vec_C=1001
        R<-object$R
        K<-object$K
        N<-object$N
        C<-object$C
        alpha<-object$optimal_model_results$alpha
        delta<-object$optimal_model_results$delta
        tau<-object$optimal_model_results$tau

        density_vec<-vector()
        cluster_vec<-vector()
        sample_vec<-vector()
        beta_vec<-vector()
        vec_x<-seq(0.001, 0.999, length=vec_C)
        #sample_name<-c("Sample A","Sample B")

        for(i in 1:R)
        {
          for(j in 1:K)
          {
            #j=1
            tmp_vec<-sapply(vec_x,function(x) {tau[j]*(stats::dbeta(x,alpha[j,i],delta[j,i]))})
            beta_vec<-c(beta_vec,vec_x)
            density_vec<-c(density_vec,tmp_vec)
            cluster_vec<-c(cluster_vec,rep(j,times=length(tmp_vec)))
            sample_vec<-c(sample_vec,rep(sample_name[i],times=length(tmp_vec)))

          }
        }

        df_new_tmp<-as.data.frame(cbind(beta_vec,density_vec,cluster_vec,sample_vec))
        df_new_tmp$sample_vec<-as.factor(df_new_tmp$sample_vec)
        df_new_tmp$cluster_vec<-as.factor(df_new_tmp$cluster_vec)
        df_new_tmp$beta_vec<-as.numeric(df_new_tmp$beta_vec)
        df_new_tmp$density_vec<-as.numeric(df_new_tmp$density_vec)
        color_length<-R
        colours<-scales::seq_gradient_pal(low="#FFC20A",high="#0C7BDC",space = "Lab")(1:color_length/color_length)

        plot_graph<-ggplot2::ggplot(df_new_tmp,ggplot2::aes(x=beta_vec,y=density_vec,color=sample_vec))+
          ggplot2::geom_line()+
          ggplot2::scale_color_manual(values=colours)+
          ggplot2::facet_wrap(~cluster_vec,scales = scale_param
          )+ ggplot2::labs(color="DNA sample", x="Beta value", y="Density")+
          ggplot2::ggtitle(txt)
        f_labels<-data.frame(cluster_vec=as.factor(seq(1,K,by=1)),label=as.vector(round((object$optimal_model_results$cluster_size/object$C),3)))
        plot_graph<-plot_graph+
          ggplot2::geom_text(data = f_labels, ggplot2::aes(x = 0.2, y = 0.1,label = label,color=NA),show.legend = F,fontface="bold" )

      }

    }

  }
  if(threshold==TRUE)
  {
    if(!is.null(plot_graph)){
      if(object$optimal_model == "K..")
      {
        th_plot<-object$optimal_model_results$thresholds
        ano_th<-object$optimal_model_results$thresholds-.02
      }else{
        th_plot<-object$optimal_model_results$thresholds$threholds[,pn]
        ano_th<-object$optimal_model_results$thresholds$threholds[,pn]-.02
      }
      num=sort(p.text$y)
      ano_y=num[length(num)]-0.1
      plot_graph<-plot_graph+ggplot2::geom_vline(xintercept = th_plot,linetype="dotted")+ggplot2::annotate("text",x=ano_th,y=ano_y,label=th_plot,angle=90)
    }
  }

  if(what == "uncertainty")
  {
    labels<-c(max_uc="Maximum uncertainty")
    #unc_df<-cbind(object$optimal_model_results$uncertainty,object$optimal_model_results$data[,"mem_final"])
    unc_df<-cbind(object$optimal_model_results$uncertainty,object$optimal_model_results$classification)
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
      )+ggplot2::labs(color="Cluster number")+
      ggplot2::xlab("Cluster number")+
      #ggplot2::theme(legend.position = "none")+
      ggplot2::ggtitle(txt)+
      #ggplot2::ggtitle("Boxplot for uncertainties in clustering solution")+
      ggplot2::coord_cartesian(ylim = c(0, 1))+
      ggplot2::geom_hline(ggplot2::aes(yintercept=max_unc, linetype="Maximum uncertainty"),color="black")+
      ggplot2::scale_linetype_manual(name="", values = 2,guide = ggplot2::guide_legend(override.aes = list(color = "black")))



  }
  if(what == "information criterion")
  {
    if(length(object$ic_output)>1)
    {
      model_name<-as.vector(object$function_call$model_names[2:length(object$function_call$model_names)])
      model_name_wo_dot<-gsub("[[:punct:]]", "", model_name)
      for(i in 1:length(model_name_wo_dot)){
        if(model_name_wo_dot[i]=="K")
        {
          model_name_wo_dot[i]="K.."
        }else if(model_name_wo_dot[i]=="KN")
        {
          model_name_wo_dot[i]="KN."
        }else if(model_name_wo_dot[i]=="KR")
        {
          model_name_wo_dot[i]="K.R"
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


  if(!is.null(plot_graph))
  {
    if(plot_type=="ggplot")
      plot_graph
    else
      plotly::ggplotly(plot_graph)
  }


}
