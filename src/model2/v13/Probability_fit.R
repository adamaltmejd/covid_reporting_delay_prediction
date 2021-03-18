
probability_analysis <- function(model_death_dt ,days.reported=5, lag = 20){

    report_cleaned <- report_clean(model_death_dt$detected,model_death_dt$dates)
    new_cases <- newCases(report_cleaned)
    rownames(new_cases) <- as.character(model_death_dt$dates)
    colnames(new_cases) <- as.character(model_death_dt$dates)
    diag(new_cases) <- NA
    N <- apply(new_cases,1,sum, na.rm=T)
    Prob <- matrix(0,nrow=dim(new_cases)[1]-lag,ncol=30)
    Analysis_data <- c()
    coeff <- c()
    dates <- c()
    for(i in (days.reported+1):(dim(new_cases)[1]-lag)){
        cases <- sum(new_cases[(i-days.reported),1:i], na.rm=T)
        Total <- sum(new_cases[i-days.reported,], na.rm=T)
        if(is.na(cases)==F & Total >0){


            Analysis_data <- rbind(Analysis_data, c(cases, Total))
            Analysis_glm  <- data.frame(Y = Analysis_data[,1],
                                        N = Analysis_data[,2]-Analysis_data[,1])
            fit <- glm(cbind(Y,N)~1, data=Analysis_glm, family=binomial)
            coeff <- c(coeff,fit$coefficients)
            dates <- c(dates,colnames(new_cases)[i])
        }
    }

    plot.data <- data.frame(date = dates, prob = 1/(1+exp(-coeff)))
    ggfig <- ggplot(data = plot.data, aes(y = prob, x = date)) +
        geom_point()
    return(list(fit = plot.data, fig = ggfig, new_cases= new_cases))
}

graphics.off()
model_death_dt <- readRDS(file.path("data", "processed", "processed_data.rds"))
days.reported <- 10

res <- probability_analysis(model_death_dt,days.reported=days.reported)
print(res$fig)


i <- dim(res$new_cases)[1]
cases <- sum(res$new_cases[(i-days.reported),1:i], na.rm=T)
size <- 0:200
prob_size = dbinom(cases,size,prob=res$fit[dim(res$fit)[1],2])
plot(size, prob_size/sum(prob_size),typ='l')
