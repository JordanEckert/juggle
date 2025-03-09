#### Simulation Study ####

# Libraries
library(utils)
devtools::install_github("JordanEckert/juggle")
library(juggle)

setwd("~/DataspellProjects/EP3CD - DELETE AFTER SIM STUDY/simulations/simulatedate_experiments")

# Dataset Generation
generate_data <- function(n, m, d, delta) {
  # Simulate data from first class
  class1_data <- matrix(runif(n * d), nrow = n)

  # Simulate data from second class
  class2_data <- matrix(runif(m * d, min = 0 + delta, max = 1 + delta), nrow = m)

  # Combine data from both classes
  data <- rbind(class1_data, class2_data)

  # Create numeric labels
  labels_numeric <- c(rep(1, n), rep(2, m))

  # Convert numeric labels to factor
  labels_factor <- factor(labels_numeric, levels = c(1, 2), labels = c("Class 1", "Class 2"))

  # Combine data and labels into a data frame
  df <- data.frame(data)
  df$Class <- labels_factor

  # Return data frame
  return(df)
}

# Inputs
dimensions <- c(2)
imbalance_ratio <- c(.1, .25, .50, .75, 1)
deltas <- c(0, .25, .50, .75, 1)

# Looping through inputs
for(dim in seq_along(dimensions)) {
  for (imbs in seq_along(imbalance_ratio)){
    for(dels in seq_along(deltas)){

      #### Parameters ####
      d <- dimensions[dim] # Number of features
      n <- 500 # Number of samples from first class
      ir <- imbalance_ratio[imbs] # Imbalance ratio
      m <- floor(n * ir) # Number of samples from second class
      delta <- deltas[dels] # Shift in the second class

      # output files
      outputfile <- paste("./results/results_dim_", d,"_delta_", sprintf("%.2f",delta), "_ir_", sprintf("%.2f",ir), ".txt",sep="")
      countfile <- paste("./results/counts_dim_", d,"_delta_",sprintf("%.2f",delta),"_ir_", sprintf("%.2f", ir),".txt",sep="")

      finaltable <- NULL
      counttable <- NULL

      # parameter sets
      e <- seq(0,1,0.1)
      tau <- seq(0,1,0.1)
      tau[1] <- .Machine$double.eps
      k <- 1:30

      # parameter counts
      inde <- rep(0,length(e))
      indtau <- rep(0,length(tau))

      indtau.acdc1 <- rep(0,length(tau))
      indtau.acdc2 <- rep(0,length(tau))
      indtau.acdc3 <- rep(0,length(tau))
      indtau.acdc4 <- rep(0,length(tau))
      indtau.acdc5 <- rep(0,length(tau))

      ## Pilot study
      for(j in 1:200){
        print(paste("Dimension=",dimensions[dim],sep=""))
        print(paste("Imbalance=",imbalance_ratio[imbs],sep=""))
        print(paste("Deltas=",deltas[dels],sep=""))
        print(paste("Pilot=",j,sep=""))

        # matrix of simulation result
        listpremc <- NULL
        listauc <- NULL

        train_data <- generate_data(n, m, d, delta)
        test_data <- generate_data(n, m, d, delta)

        ## P-CCCD
        premc <- rep(0,length(tau))

        for(z in 1:length(tau)){
          model <- juggle::pcccd(train_data[,-ncol(train_data)], as.factor(train_data$Class), tau = tau[z])
          results <- juggle::classify_pcccd(model, test_data[,-ncol(test_data)])

          premc[z] <- juggle::auc(as.numeric(results), as.numeric(test_data$Class))
        }

        tp <- which(premc==max(premc))
        indtau[tp] <- indtau[tp] + 1

        ## RW-CCCD
        premc <- rep(0,length(e))

        for(z in 1:length(e)){
          model <- juggle::rwcccd(train_data[,-ncol(train_data)], as.factor(train_data$Class))
          results <- juggle::classify_rwcccd(model, test_data[,-ncol(test_data)], e = e[z])

          premc[z] <- juggle::auc(as.numeric(results), as.numeric(test_data$Class))
        }

        tp <- which(premc==max(premc))
        inde[tp] <- inde[tp] + 1

        ## ACDC - 1
        premc <- rep(0,length(tau))

        for(z in 1:length(tau)){
          model <- juggle::acdc(train_data[,-ncol(train_data)], as.factor(train_data$Class),
                                tau = tau[z], test_data[,-ncol(test_data)], num_iter = 10, num_balls = 1)
          results <- as.numeric(as.factor(model$final_predictions))

          premc[z] <- juggle::auc(as.numeric(results), as.numeric(test_data$Class))
        }

        indtau.acdc1[tp] <- indtau.acdc1[tp] + 1

        ## ACDC - 2
        premc <- rep(0,length(tau))

        for(z in 1:length(tau)){
          model <- juggle::acdc(train_data[,-ncol(train_data)], as.factor(train_data$Class),
                                tau = tau[z], test_data[,-ncol(test_data)], num_iter = 10, num_balls = 2)
          results <- as.numeric(as.factor(model$final_predictions))

          premc[z] <- juggle::auc(as.numeric(results), as.numeric(test_data$Class))
        }

        indtau.acdc2[tp] <- indtau.acdc2[tp] + 1

        ## ACDC - 3
        premc <- rep(0,length(tau))

        for(z in 1:length(tau)){
          model <- juggle::acdc(train_data[,-ncol(train_data)], as.factor(train_data$Class),
                                tau = tau[z], test_data[,-ncol(test_data)], num_iter = 10, num_balls = 3)
          results <- as.numeric(as.factor(model$final_predictions))

          premc[z] <- juggle::auc(as.numeric(results), as.numeric(test_data$Class))
        }

        indtau.acdc3[tp] <- indtau.acdc3[tp] + 1

        ## ACDC - 4
        premc <- rep(0,length(tau))

        for(z in 1:length(tau)){
          model <- juggle::acdc(train_data[,-ncol(train_data)], as.factor(train_data$Class),
                                tau = tau[z], test_data[,-ncol(test_data)], num_iter = 10, num_balls = 4)
          results <- as.numeric(as.factor(model$final_predictions))

          premc[z] <- juggle::auc(as.numeric(results), as.numeric(test_data$Class))
        }

        indtau.acdc4[tp] <- indtau.acdc4[tp] + 1

        ## ACDC - 5
        premc <- rep(0,length(tau))

        for(z in 1:length(tau)){
          model <- juggle::acdc(train_data[,-ncol(train_data)], as.factor(train_data$Class),
                                tau = tau[z], test_data[,-ncol(test_data)], num_iter = 10, num_balls = 5)
          results <- as.numeric(as.factor(model$final_predictions))

          premc[z] <- juggle::auc(as.numeric(results), as.numeric(test_data$Class))
        }

        indtau.acdc5[tp] <- indtau.acdc5[tp] + 1
      }

        #### Writing Parameter Count Files for ACDC ####
        count1 <- c(d, delta, ir, indtau.acdc1)
        count2 <- c(d, delta, ir, indtau.acdc2)
        count3 <- c(d, delta, ir, indtau.acdc3)
        count4 <- c(d, delta, ir, indtau.acdc4)
        count5 <- c(d, delta, ir, indtau.acdc5)

        countable <- rbind(count1, count2, count3, count4, count5)
        write.table(countable, file=countfile)

        #### Finding Optimal Parameters from Pilot ####
        opt.e <-  e[which.max(inde)]
        opt.tau <-  tau[which.max(indtau)]
        opt.tau.1 <- tau[which.max(indtau.acdc1)]
        opt.tau.2 <- tau[which.max(indtau.acdc2)]
        opt.tau.3 <- tau[which.max(indtau.acdc3)]
        opt.tau.4 <- tau[which.max(indtau.acdc4)]
        opt.tau.5 <- tau[which.max(indtau.acdc5)]

        #### Monte Carlo ####
        niter <- 2500
        cccds.auc <- rep(0,2)
        ecccds.auc <- rep(0,5)

        for(j in 1:niter){

          # current iter
          print(paste("Iter=",j,sep=""))

          train_data <- generate_data(n, m, d, delta)
          test_data <- generate_data(n, m, d, delta)

          ## P-CCCD
          pmodel <- juggle::pcccd(x = train_data[,-ncol(train_data)],
                         y = as.factor(train_data[,ncol(train_data)]),
                         tau = opt.tau)
          presults <- juggle::classify_pcccd(pmodel, test_data[,-ncol(test_data)])
          cccds.auc[1] <- juggle::auc(as.numeric(presults), as.numeric(test_data$Class))


          ## RW-CCCD
          rwmodel <- juggle::rwcccd(x = train_data[,-ncol(train_data)],
                         y = as.factor(train_data[,ncol(train_data)]))
          rwresults <- juggle::classify_rwcccd(rwmodel, test_data[,-ncol(test_data)], e = opt.e)
          cccds.auc[2] <- juggle::auc(as.numeric(rwresults), as.numeric(test_data$Class))

          ## ACDC - 1
          a1model <- juggle::acdc(train_data[,-ncol(train_data)], as.factor(train_data$Class),
                                tau = opt.tau.1, test_data[,-ncol(test_data)], num_iter = 10, num_balls = 1)
          a1results <- as.numeric(as.factor(a1model$final_predictions))

          ecccds.auc[1] <- juggle::auc(as.numeric(a1results), as.numeric(test_data$Class))

          ## ACDC - 2
          a2model <- juggle::acdc(train_data[,-ncol(train_data)], as.factor(train_data$Class),
                                tau = opt.tau.2, test_data[,-ncol(test_data)], num_iter = 10, num_balls = 2)
          a2results <- as.numeric(as.factor(a2model$final_predictions))

          ecccds.auc[2] <- juggle::auc(as.numeric(a2results), as.numeric(test_data$Class))

          ## ACDC - 3
          a3model <- juggle::acdc(train_data[,-ncol(train_data)], as.factor(train_data$Class),
                                tau = opt.tau.3, test_data[,-ncol(test_data)], num_iter = 10, num_balls = 3)
          a3results <- as.numeric(as.factor(a3model$final_predictions))

          ecccds.auc[3] <- juggle::auc(as.numeric(a3results), as.numeric(test_data$Class))

          ## ACDC - 4
          a4model <- juggle::acdc(train_data[,-ncol(train_data)], as.factor(train_data$Class),
                                tau = opt.tau.4, test_data[,-ncol(test_data)], num_iter = 10, num_balls = 4)
          a4results <- as.numeric(as.factor(a4model$final_predictions))

          ecccds.auc[4] <- juggle::auc(as.numeric(a4results), as.numeric(test_data$Class))

          ## ACDC - 5
          a5model <- juggle::acdc(train_data[,-ncol(train_data)], as.factor(train_data$Class),
                                tau = opt.tau.5, test_data[,-ncol(test_data)], num_iter = 10, num_balls = 5)

          a5results <- as.numeric(as.factor(a5model$final_predictions))

          ecccds.auc[5] <- juggle::auc(as.numeric(a5results), as.numeric(test_data$Class))

          max.auc <- c(cccds.auc,ecccds.auc)
          listauc <- rbind(listauc,max.auc)

          if(j==1) print(max.auc)
          else{
            ave.auc <- apply(listauc,2,mean)
            sd.auc <- apply(listauc,2,sd)/(j)

            print(rbind(ave.auc,sd.auc))

            if(all(sd.auc < 0.001)) break
          }
        }
        result <- c(d, delta, ir, ave.auc)
        finaltable <- rbind(finaltable,result)

        write.table(finaltable, file=outputfile)
    }
  }
}
