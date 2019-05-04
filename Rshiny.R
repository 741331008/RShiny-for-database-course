source("Rfunc.R")

# Define ui
ui <- fluidPage(
    
    # Application title
    titlePanel("Data Visualization"),

    # Siderbar layout
    sidebarLayout(
        # show the input and design of sidebar
        sidebarPanel(
            selectInput("type","Cancer Type:",choices = c("Acute Myeloid Leukemia","Non Small Cell Lung Cancer","Breast Cancer")),
            selectInput("variable","Grouping Variable:",choices = c("Samples","Genes","Mutation Types")),
            selectInput("dist","Distance between samples:",choices = c("euclidean","manhattan","maximum","binary"))
        ),

        # show the output and different tables.
        mainPanel(
            tabsetPanel(type = "tabs",
                        tabPanel("Statistics",dataTableOutput("Samples")),
                        tabPanel("Signatures",dataTableOutput("Signatures")),
                        tabPanel("Actionable_Table",dataTableOutput("Guidelines")),
                        tabPanel("Heatmap",d3heatmapOutput("heatmap",height="1000px")),
                        tabPanel("Violin Plot",plotOutput("violinPlot",height = "600px"))
            )
        )
    )
)

shinyServer <- function(input,output,session){
    
    myConn = dbConnect(SQLite(),
                       "C:\\Users\\Super Kang\\Documents\\class\\2019 Spring\\Database Management\\final project\\pharmcogenome.db")
    output$Guidelines <- renderDataTable({dbGetQuery(myConn,
                                    "SELECT * FROM allActionableVariants")})
    
    output$Samples <- renderDataTable({database_cal(input$type,input$variable)})
    output$Signatures <- renderDataTable({survival_relatedMutations(input$type)})
    output$heatmap <- renderD3heatmap(d3heatmap(generate_heatmap(input$type),
                                                colors = "Blues",
                                                distfun = function(x) dist(x,method = input$dist)))

    output$violinPlot <- renderPlot({
        if (input$type=="Non Small Cell Lung Cancer")
        {
            tbl <- violin_forSignatures(input$type,5)
        p1 <- ggplot(tbl,aes(x=tbl[,3],y=Survival,fill=tbl[,3])) +geom_violin(alpha=0.6) + theme(legend.position = "none") + xlab(colnames(tbl)[3]) + ylab("ProgressFree Survival (Months)")
        p2 <- ggplot(tbl,aes(x=tbl[,4],y=Survival,fill=tbl[,4])) +geom_violin(alpha=0.6) + theme(legend.position = "none") + xlab(colnames(tbl)[4]) + ylab("ProgressFree Survival (Months)")
        p3 <- ggplot(tbl,aes(x=tbl[,5],y=Survival,fill=tbl[,5])) +geom_violin(alpha=0.6) + theme(legend.position = "none") + xlab(colnames(tbl)[5]) + ylab("ProgressFree Survival (Months)")
        grid.arrange(p1,p2,p3,ncol=2)}
        if (input$type=="Breast Cancer")
        {
            tbl <- violin_forSignatures(input$type)
            p1 <- ggplot(tbl,aes(x=tbl[,3],y=Survival,fill=tbl[,3])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[3]) + ylab("Survival (Months)")
            p2 <- ggplot(tbl,aes(x=tbl[,4],y=Survival,fill=tbl[,4])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[4]) + ylab("Survival (Months)")
            p3 <- ggplot(tbl,aes(x=tbl[,5],y=Survival,fill=tbl[,5])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[5]) + ylab("Survival (Months)")
            p4 <- ggplot(tbl,aes(x=tbl[,6],y=Survival,fill=tbl[,6])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[6]) + ylab("Survival (Months)")
            p5 <- ggplot(tbl,aes(x=tbl[,7],y=Survival,fill=tbl[,7])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[7]) + ylab("Survival (Months)")
            p6 <- ggplot(tbl,aes(x=tbl[,8],y=Survival,fill=tbl[,8])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[8]) + ylab("Survival (Months)")
            p7 <- ggplot(tbl,aes(x=tbl[,9],y=Survival,fill=tbl[,9])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[9]) + ylab("Survival (Months)")
            p8 <- ggplot(tbl,aes(x=tbl[,10],y=Survival,fill=tbl[,10])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[10]) + ylab("Survival (Months)")
            p9 <- ggplot(tbl,aes(x=tbl[,11],y=Survival,fill=tbl[,11])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[11]) + ylab("Survival (Months)")
            p10 <- ggplot(tbl,aes(x=tbl[,12],y=Survival,fill=tbl[,12])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[12]) + ylab("Survival (Months)")
            p11 <- ggplot(tbl,aes(x=tbl[,13],y=Survival,fill=tbl[,13])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[13]) + ylab("Survival (Months)")
            p12 <- ggplot(tbl,aes(x=tbl[,14],y=Survival,fill=tbl[,14])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[14]) + ylab("Survival (Months)")
            p13 <- ggplot(tbl,aes(x=tbl[,15],y=Survival,fill=tbl[,15])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[15]) + ylab("Survival (Months)")
            p14 <- ggplot(tbl,aes(x=tbl[,16],y=Survival,fill=tbl[,16])) +geom_violin(alpha=0.6) + theme(legend.position = "none") +xlab(colnames(tbl)[16]) + ylab("Survival (Months)")
            grid.arrange(p1,p2,p3,p4,p5,p6,p7,p8,p9,p10,p11,p12,p13,p14,ncol=3)
        }
    })
    
    
}

shinyApp(ui=ui,server=shinyServer)