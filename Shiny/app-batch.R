

library(MASS)
library(shiny)
library(DT)
library(shinythemes)
library(waiter)

# library(rsconnect)
# setRepositories()
# options(repos = BiocManager::repositories())
# deployApp("/home/mquinodo/SYNO/MutLand/Shiny-MutScore",appName="mutscore-batch")

# deployApp("~/Desktop/Shiny-MutScore",appName="mutscore-batch")
# shiny::runApp("/home/mquinodo/SYNO/MutLand/Shiny-MutScore")

load(file="genelist.RData")
load(file="PLP-BLB.RData")

ui <- navbarPage(id = "start",
  windowTitle="MutScore",
  
  #waiter_show_on_load(spin_throbber()),
  
  theme = shinytheme("cerulean"),
  title=div(img(src="IOB2.png", height="30px", width="110px"), style = "font-size:30px",strong("MutScore")),
  tabPanel(value="end",title = "Single variant",
    fluidRow(
      sidebarPanel(
        selectizeInput('chr', "Chromosome", as.list(genelist),selected="1"),
        numericInput('position', 'Position', 69091,min = 1, max = 10000000000),
        selectizeInput('ref', 'Reference', as.list(c("A","C","T","G")),selected="A"),
        selectizeInput('alt', 'Alternative', as.list(c("A","C","T","G")),selected="G"),
        downloadButton('download1', 'Download data as tsv'),
        selectizeInput('build', 'Genome build', as.list(c("","hg19","hg38")),selected=""),
        width=2
      ),
      mainPanel(DT::dataTableOutput('table'),br(),plotOutput('plot1',height="550px",width="950px"),width=8)
    )
  ),
  
  tabPanel(value="start",title = "Genomic region",
    #use_waiter(), 
           fluidRow(
             sidebarPanel(
               selectizeInput('chr2', "Chromosome", as.list(genelist),selected="1"),
               numericInput('begin', 'Begin', 69091,min = 1, max = 10000000000),
               numericInput('end', 'End', 69094,min = 1, max = 10000000000),
               downloadButton('download2', 'Download data as tsv'),
               width=2
             ),
             column(9,offset=0,
                    DT::dataTableOutput('table2')
             )
           )
  ),
  
  tabPanel(title = "List of variants (batch query)",
           use_waiter(),
           use_steward(),
           fluidRow(
             sidebarPanel(
               textAreaInput('list', "List of variants\ntest", value = "chr1:897020G>C\n1:897012A>T",rows=10, placeholder = NULL, resize = NULL),
               actionButton("button", "Refresh table"),
               downloadButton('download3', 'Download data as tsv'),
               width=2
             ),
             column(8,offset=0,
                    DT::dataTableOutput('table3')
             )
           )
  ),
  
  tabPanel(title = "About - contact",
           titlePanel(h3("About")),
           p("Large-scale sequencing efforts of genomes from patients often lead to the identification of DNA variants of unknown significance (VUS), i.e. of genomic variations that cannot be immediately recognized as pathogenic mutations or benign DNA changes. Despite many in silico predictors exist, none of them takes advantage from the wealth of information contained in the positional clustering of mutations already detected in disease-associated genes."),
             p("We have therefore developed MutScore, a new pathogenicity predictor that integrates unsupervised features of single DNA variants with information derived from such clusters. The predictive model for MutScore was trained with a random forest approach on medically-relevant mutations and subsequently tested against various genomic databases for both hereditary conditions and cancer (ClinVar, HGMD, and DoCM), achieving a very high performance."),
             p("The use of MutScore on 840 genes from the ClinVar database also allowed the detection of significant clusters of disease-associated and of benign variants in 505 and 345 of them, respectively, revealing protein domains with diverging functional importance. In addition, an open-access web-based application, MutLand, was developed to provide a comprehensive graphical landscape of all known medically-relevant and clearly benign DNA variants for individual genes, as a help in appraising new VUS identified in clinical testing."),
             p("Altogether, our work reveals the presence of widespread clustering of missense variants associated with normal and clinical phenotypes and that this information can be systematically used to improve and to understand pathogenicity at the molecular level."),
           br(),
           titlePanel(h3("Contact")),
           p("Mathieu Quinodoz"),
           div("mathieu.quinodoz[at]iob.ch", style = "color:blue"),
           br(),
           p("Institute of Molecular and Clinical Ophthalmology Basel (IOB)"),
           p("Mittlere Strasse 91, 4031 Basel, Switzerland"),
           a("IOB website", href="https://www.iob.ch/"),
           br(),
           img(src = "IOB2.png", height="90px", width="330px")
  ),
  
  tabPanel(title = "Download",
           titlePanel(h3("MutScore")),
           br(),
           p("MutScore and MutLand are freely available for non-commercial use. For other uses, please contact IOB."),
           br(),
           p("MutScore values for all missense variants for hg19 genome build can be downloaded on our github page:"),
           a("github.com/mquinodo/MutScore", href="https://github.com/mquinodo/MutScore", target="_blank")
  ),
  
  tabPanel(title = "How to cite",
           titlePanel(h3("Manuscript")),
           br(),
           p("A manuscript about MutScore has just been submitted."),
  ),
  
  tabPanel(title = div(strong("MutScore"),style = "font-size:20px"),
           # div(img(src="IOB2.png", height="30px", width="110px"), style = "font-size:30px",strong("MutScore"))
           titlePanel(h3("MutScore")),
           br(),
           p("For MutScore, please visit this webpage:"),
           a("iob-genetic.shinyapps.io/mutscore", href="https://iob-genetic.shinyapps.io/mutscore", target="_blank")
  )
)

server <- function(input, output, session) {
  
  
  load(file="data12-mutscore-test.RData")
  
  v <- reactiveValues()
  
  showModal(modalDialog(
    title = "Choose reference genome:",
    actionButton("button19", "GRCh37/hg19"),actionButton("button38", "GRCh38/hg38")
  ))
  
  v$tem=""
  
  observeEvent(input$button19, {
    removeModal()
    waiter_show(html = span("Please wait while loading data for hg19 (less than 1 minute)"))
    #Sys.sleep(10)
    v$mutscore2=readRDS(file="data12-mutscore.RDS")
    colnames(v$mutscore2)=c("Chromosome","Position","Reference","Alternative","MutScore")
    waiter_hide()
    updateTextInput(session,"build",value="hg19")
    v$tem="hg19"
  })
  
  observeEvent(input$button38, {
    removeModal()
    waiter_show(html = span("Please wait while loading data for hg38 (less than 1 minute)"))
    print("a")
    v$mutscore2=readRDS(file="data12-mutscore-hg38.RDS")
    print("b")
    colnames(v$mutscore2)=c("Chromosome","Position","Reference","Alternative","MutScore")
    print("c")
    waiter_hide()
    updateTextInput(session,"build",value="hg38")
    v$tem="hg38"
  })
  
  observeEvent(input$build, {
    if(input$build=="hg19" & v$tem!="hg19"){
      waiter_show(html = span("Please wait while loading data for hg19 (less than 1 minute)"))
      v$mutscore2=readRDS(file="data12-mutscore.RDS")
      colnames(v$mutscore2)=c("Chromosome","Position","Reference","Alternative","MutScore")
      waiter_hide()
    }
    if(input$build=="hg38" & v$tem!="hg38"){
      waiter_show(html = span("Please wait while loading data for hg38 (less than 1 minute)"))
      v$mutscore2=readRDS(file="data12-mutscore-hg38.RDS")
      colnames(v$mutscore2)=c("Chromosome","Position","Reference","Alternative","MutScore")
      waiter_hide()
    }
  })
  
  #w1 <- Waiter$new(spin_throbber(),id = "table")
  w2 <- Waiter$new(spin_throbber(),id = "table2")

  output$table <- DT::renderDataTable({
    #w1$show()
    #on.exit({
    #  w1$hide()
    #})
    a=which(v$mutscore2[,1]==input$chr)
    b=which(v$mutscore2[a,2]==input$position)
    c=which(v$mutscore2[a[b],3]==input$ref & v$mutscore2[a[b],4]==input$alt)
    DT::datatable(v$mutscore2[a[b[c]],], options = list(paging=FALSE,searching=FALSE,columnDefs = list(list(className = 'dt-center', targets=0:4))),rownames=FALSE)
  })
  
  output$table2 <- DT::renderDataTable({
    w2$show()
    on.exit({
      w2$hide()
    })
    a=which(v$mutscore2[,1]==input$chr2)
    b=which(v$mutscore2[a,2]>=input$begin & v$mutscore2[a,2]<=input$end)
    DT::datatable(v$mutscore2[a[b],], options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 10,columnDefs = list(list(className = 'dt-center', targets=0:4))),rownames=FALSE)
  })
  
  w <- Waiter$new(spin_throbber(),id = "table3")
  
 out <- reactive({
    input$button
    w$show()
    on.exit({
      w$hide()
    })
    list2=unique(strsplit(gsub("chr","",isolate(input$list)),'\\n')[[1]])
    list3=""
    
    var1=""
    var2=""
    var3=""
    var4=""
    
    for (i in 1:length(list2)){
      temp1=strsplit(list2[i],':')[[1]][1]
      if(is.na(temp1)==F){
      temp2=gsub("[^0-9.-]", "",strsplit(list2[i],':')[[1]][2])
      if(is.na(temp2)==F){
      temp3=substr(gsub("[0-9.-]", "",strsplit(list2[i],':')[[1]][2]),1,1)
      if(is.na(temp3)==F){
      temp4=substr(gsub("[0-9.-]", "",strsplit(list2[i],':')[[1]][2]),3,3)
      if(is.na(temp4)==F){
      if(nchar(gsub("[0-9>]", "",strsplit(list2[i],':')[[1]][2]))==2){
        var1=c(var1,temp1)
        var2=c(var2,temp2)
        var3=c(var3,temp3)
        var4=c(var4,temp4)
      }}}}}
    }
    var1=var1[-1]
    var2=var2[-1]
    var3=var3[-1]
    var4=var4[-1]
    
    if (var2[1]==897020 & var2[2]==897012 & length(var2)==2){
      a=c(14109,14110,14111,14131,14132)
    } else {
      a=which(is.element(v$mutscore2[,2],var2))
    }
    
    mutscore3=v$mutscore2[a,]
    list3=""
    non=c("","","","","")
    for (i in 1:length(var3)){
      c=which(mutscore3[,1]==var1[i] & mutscore3[,2]==var2[i] & mutscore3[,3]==var3[i] & mutscore3[,4]==var4[i])
      if(length(c)>0){list3=c(list3,c)} else {non=rbind(non,c(var1[i],var2[i],var3[i],var4[i],"NA"))}
    }
    list3=as.numeric(list3[-1])
    if(length(non)>5){non=non[-1,]}
    
    out=mutscore3[list3,]
    if(length(non)>1){out=rbind(out,non)}
    out=out[order(as.numeric(out[,1]),as.numeric(out[,2])),]

    if(length(list3)>0){
      
    } else {
      out=c()
    }
    if(length(which(out[,5]==""))>0){out=out[-which(out[,5]==""),]}
    out
  })

  output$table3 <- DT::renderDataTable({
    input$tabs
    DT::datatable(mutscore[1:2,], options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 10,columnDefs = list(list(className = 'dt-center', targets=0:4))),rownames=FALSE)
  })
 
  output$table3 <- DT::renderDataTable({
    DT::datatable(out(), options = list(lengthMenu = c(10, 25, 50, 100), pageLength = 10,columnDefs = list(list(className = 'dt-center', targets=0:4))),rownames=FALSE)
  })

  output$download1 <- downloadHandler(
    filename = function() { paste("MutScore-",input$build,"_",input$chr,"-",input$position,input$ref,">",input$alt,'.tsv', sep='') },
    content = function(file) {
      a=which(v$mutscore2[,1]==input$chr)
      b=which(v$mutscore2[a,2]==input$position)
      c=which(v$mutscore2[a[b],3]==input$ref & v$mutscore2[a[b],4]==input$alt)
      write.table(v$mutscore2[a[b[c]],], file, row.names = FALSE,sep="\t",quote=FALSE)
    }
  )
  
  output$download2 <- downloadHandler(
    filename = function() { paste("MutScore-",input$build,"_",input$chr,"-",input$begin,"-",input$end,'.tsv', sep='') },
    content = function(file) {
      a=which(v$mutscore2[,1]==input$chr2)
      b=which(v$mutscore2[a,2]>=input$begin & v$mutscore2[a,2]<=input$end)
      write.table(v$mutscore2[a[b],], file, row.names = FALSE,sep="\t",quote=FALSE)
    }
  )

  output$download3 <- downloadHandler(
    filename = function() { paste("MutScore-",input$build,"_download_list.tsv",sep="") },
    content = function(file) {
      write.table(out(), file, row.names = FALSE,sep="\t",quote=FALSE)
    }
  )
  
  output$plot1 <- renderPlot({
    a=which(v$mutscore2[,1]==input$chr)
    b=which(v$mutscore2[a,2]==input$position)
    c=which(v$mutscore2[a[b],3]==input$ref & v$mutscore2[a[b],4]==input$alt)
    score=v$mutscore2[a[b[c]],5]
    par(mfrow=c(2,1))
    hist(PLP,n=50,xlab="MutScore",cex.main=1.5,cex.lab=1.2,col='red',main="MutScore distribution for PLP missense variants from ClinVar",yaxt='n')
    axis(2,at=c(0,3000,6000,9000,12000),labels=c(0,3000,6000,9000,12000),las=2,cex.axis=0.6)
    abline(v=score,lwd=2,col='green')
    legend(0.4, 12500, legend=c("Selected variant"),col=c("green"), lty=1, cex=0.8)
    
    n=50
    a=hist(PLP,n=n,plot=F)
    b=hist(BLB,n=n,plot=F)
    m=which.min(abs(a$mids-score))[1]
    proba=round(a$density[m]/(a$density[m]+b$density[m])*100,digits=1)
    text(0.5,8000,paste("",proba,"% of variants with similar scores are PLPs",sep=""))
    proba2=round(length(which(PLP>score))/length(PLP)*100,digits=1)
    text(0.5,6000,paste("",proba2,"% of PLPs have a score higher than ",score,sep=""))
    
    hist(BLB,n=50,xlab="MutScore",cex.main=1.5,cex.lab=1.2,col='blue',main="MutScore distribution for BLB missense variants from ClinVar",yaxt='n')
    axis(2,at=c(0,1000,2000,3000),labels=c(0,1000,2000,3000),las=2,cex.axis=0.6)
    abline(v=score,lwd=2,col='green')
    legend(0.4, 3800, legend=c("Selected variant"),col=c("green"), lty=1, cex=0.8)
    text(0.5,2500,paste("",100-proba,"% of variants with similar scores are BLBs",sep=""))
    proba2=round(length(which(BLB<score))/length(BLB)*100,digits=1)
    text(0.5,2000,paste("",proba2,"% of BLBs have a score lower than ",score,sep=""))

  })
  
}

shinyApp(ui = ui, server = server)

