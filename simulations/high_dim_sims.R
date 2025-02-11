# Libraries
library(juggle)

# Index
dimensions <- c(2500)
imbalance_ratio <- c(.1, .25, .50, .75, 1)
deltas <- c(0, .25, .50, .75, 1)

# Simulation to get Optimal Parameters
for(dim in seq_along(dimensions)) {
  for (imbs in seq_along(imbalance_ratio)){
    for(dels in seq_along(deltas)){

      # Iteration Counts
      print(paste("Dimension=",dimensions[dim],sep=""))
      print(paste("Imbalance=",imbalance_ratio[imbs],sep=""))
      print(paste("Deltas=",deltas[dels],sep=""))

      d <- dimensions[dim] # Number of features
      n <- 100 # Number of samples from first class
      ir <- imbalance_ratio[imbs] # Imbalance ratio
      m <- floor(n * ir) # Number of samples from second class
      delta <- deltas[dels] # Shift in the second class

      outputfile <- paste("./simulations/results_dim", d,"_delta_", sprintf("%.2f",delta), "_ir_", sprintf("%.2f",ir), ".txt",sep="")
      countfile <- paste("./simulations/counts_dim_", d,"_delta_",sprintf("%.2f",delta),"_ir_", sprintf("%.2f", ir),".txt",sep="")

      finaltable <- NULL
      counttable <- NULL

      # parameter sets
      e <- seq(0,1,0.1)
      tau <- seq(0,1,0.1)
      tau[1] <- .Machine$double.eps
      gamma <- seq(0.1,4,0.1)
      k <- 1:30

      # parameter counts
      indtau <- rep(0,length(tau))

      indtau.e1 <- rep(0,length(tau))
      indtau.e2 <- rep(0,length(tau))
      indtau.e3 <- rep(0,length(tau))
      indtau.e4 <- rep(0,length(tau))
      indtau.e5 <- rep(0,length(tau))

      indk.e1 <- rep(0,length(k))
      indk.e2 <- rep(0,length(k))
      indk.e3 <- rep(0,length(k))
      indk.e4 <- rep(0,length(k))
      indk.e5 <- rep(0,length(k))

      for(j in 1:200){
        # current iter
        print(paste("Pilot=",j,sep=""))

        # matrix of simulation result
        listpremc <- NULL
        listauc <- NULL

        train_data <- generate_unit_square(n, m, d, delta)
        test_data <- generate_unit_square(n, m, d, delta)

        # P-CCCD
        premc <- rep(0,length(tau))

        for(i in 1:length(tau)) {
         model <- pcccd(x = train_data[,-ncol(train_data)], y = train_data$Class, tau = tau[i])
         result <- classify_pcccd(model, test_data[,-ncol(test_data)])

         premc[i] <- auc(as.numeric(as.factor(result)), as.numeric(as.factor(test_data$Class)))
        }

        tp <- which(premc==max(premc))
        indtau[tp] <- indtau[tp] + 1

        # ACDC k = 1
        premc <- rep(0, length(tau))
        premc2 <- rep(0, length(k))

        for(i in 1:length(tau)){
          for (j in 1:length(k)){
            model <- acdc(x = train_data[,-ncol(train_data)], y = train_data$Class,
                          test_data = test_data[,-ncol(test_data)], tau = tau[i], num_iter = k[j],
                          num_balls = 1)
            result <- as.factor(model$final_predictions)

            premc2[j] <- auc(as.numeric(as.factor(result)), as.numeric(as.factor(test_data$Class)))
          }
          premc[i] <- max(premc2)
        }

        tp <- which(premc==max(premc))
        tp2 <- which(premc2==max(premc2))
        indtau.e1[tp] <- indtau.e1[tp] + 1
        indk.e1[tp2] <- indk.e1[tp2] + 1

        # ACDC k = 2
        premc <- rep(0, length(tau))
        premc2 <- rep(0, length(k))

        for(i in 1:length(tau)){
          for (j in 1:length(k)){
            model <- acdc(x = train_data[,-ncol(train_data)], y = train_data$Class,
                          test_data = test_data[,-ncol(test_data)], tau = tau[i], num_iter = k[j],
                          num_balls = 2)
            result <- as.factor(model$final_predictions)

            premc2[j] <- auc(as.numeric(as.factor(result)), as.numeric(as.factor(test_data$Class)))
          }
          premc[i] <- max(premc2)
        }

        tp <- which(premc==max(premc))
        tp2 <- which(premc2==max(premc2))
        indtau.e2[tp] <- indtau.e2[tp] + 1
        indk.e2[tp2] <- indk.e2[tp2] + 1

        # ACDC k = 3
        premc <- rep(0, length(tau))
        premc2 <- rep(0, length(k))

        for(i in 1:length(tau)){
          for (j in 1:length(k)){
            model <- acdc(x = train_data[,-ncol(train_data)], y = train_data$Class,
                          test_data = test_data[,-ncol(test_data)], tau = tau[i], num_iter = k[j],
                          num_balls = 3)
            result <- as.factor(model$final_predictions)

            premc2[j] <- auc(as.numeric(as.factor(result)), as.numeric(as.factor(test_data$Class)))
          }
          premc[i] <- max(premc2)
        }

        tp <- which(premc==max(premc))
        tp2 <- which(premc2==max(premc2))
        indtau.e3[tp] <- indtau.e3[tp] + 1
        indk.e3[tp2] <- indk.e3[tp2] + 1

        # ACDC k = 4
        premc <- rep(0, length(tau))
        premc2 <- rep(0, length(k))

        for(i in 1:length(tau)){
          for (j in 1:length(k)){
            model <- acdc(x = train_data[,-ncol(train_data)], y = train_data$Class,
                          test_data = test_data[,-ncol(test_data)], tau = tau[i], num_iter = k[j],
                          num_balls = 4)
            result <- as.factor(model$final_predictions)

            premc2[j] <- auc(as.numeric(as.factor(result)), as.numeric(as.factor(test_data$Class)))
          }
          premc[i] <- max(premc2)
        }

        tp <- which(premc==max(premc))
        tp2 <- which(premc2==max(premc2))
        indtau.e4[tp] <- indtau.e4[tp] + 1
        indk.e4[tp2] <- indk.e4[tp2] + 1

        # ACDC k = 5
        premc <- rep(0, length(tau))
        premc2 <- rep(0, length(k))

        for(i in 1:length(tau)){
          for (j in 1:length(k)){
            model <- acdc(x = train_data[,-ncol(train_data)], y = train_data$Class,
                          test_data = test_data[,-ncol(test_data)], tau = tau[i], num_iter = k[j],
                          num_balls = 5)
            result <- as.factor(model$final_predictions)

            premc2[j] <- auc(as.numeric(as.factor(result)), as.numeric(as.factor(test_data$Class)))
          }
          premc[i] <- max(premc2)
        }

        tp <- which(premc==max(premc))
        tp2 <- which(premc2==max(premc2))
        indtau.e5[tp] <- indtau.e5[tp] + 1
        indk.e5[tp2] <- indk.e5[tp2] + 1

        #### Writing Parameter Count Files for ACDC ####
        count1 <- c(d, delta, ir, indtau.e1, indk.e1)
        count2 <- c(d, delta, ir, indtau.e2, indk.e2)
        count3 <- c(d, delta, ir, indtau.e3, indk.e3)
        count4 <- c(d, delta, ir, indtau.e4, indk.e4)
        count5 <- c(d, delta, ir, indtau.e5, indk.e5)

        countable <- rbind(count1, count2, count3, count4, count5)
        colnames(countable) <- c("Dimension", "Delta", "Imbalance", "Tau1", "Tau2", "Tau3", "Tau4", "Tau5",
                                 "Tau6", "Tau7", "Tau8", "Tau9", "Tau10", "Tau11", "K1", "K2", "K3", "K4", "K5",
                                 "K6", "K7", "K8", "K9", "K10", "K11", "K12", "K13", "K14", "K15", "K16", "K17",
                                 "K18", "K19", "K20", "K21", "K22", "K23", "K24", "K25", "K26", "K27", "K28", "K29", "K30")
        write.table(countable, file=countfile)

        #### Finding Optimal Parameters from Pilot ####
        opt.tau <-  tau[which.max(indtau)]
        opt.tau.1 <- tau[which.max(indtau.e1)]
        opt.tau.2 <- tau[which.max(indtau.e2)]
        opt.tau.3 <- tau[which.max(indtau.e3)]
        opt.tau.4 <- tau[which.max(indtau.e4)]
        opt.tau.5 <- tau[which.max(indtau.e5)]

        opt.k.e1 <-  k[which.max(indk.e1)]
        opt.k.e2 <-  k[which.max(indk.e2)]
        opt.k.e3 <-  k[which.max(indk.e3)]
        opt.k.e4 <-  k[which.max(indk.e4)]
        opt.k.e5 <-  k[which.max(indk.e5)]

        #### Monte Carlo ####
        niter <- 2500
        cccd.auc <- rep(0,1)
        ecccds.auc <- rep(0,5)

        for(j in 1:niter){
          # current iter
          print(paste("Iter=",j,sep=""))

          train_data <- generate_unit_square(n, m, d, delta)
          test_data <- generate_unit_square(n, m, d, delta)

          # P-CCCD
          pmodel <- pcccd(x = train_data[,-ncol(train_data)], y = train_data$Class, tau = opt.tau)
          result <- classify_pcccd(pmodel, test_data[,-ncol(test_data)])
          cccd.auc[1] <- auc(as.numeric(as.factor(result)), as.numeric(as.factor(test_data$Class)))

          # ACDC k = 1
          model <- acdc(x = train_data[,-ncol(train_data)], y = train_data$Class,
                        test_data = test_data[,-ncol(test_data)], tau = opt.tau.1, num_iter = opt.k.e1,
                        num_balls = 1)
          result <- as.factor(model$final_predictions)
          ecccds.auc[1] <- auc(as.numeric(as.factor(result)), as.numeric(as.factor(test_data$Class)))

          # ACDC k = 2
          model2 <- acdc(x = train_data[,-ncol(train_data)], y = train_data$Class,
                        test_data = test_data[,-ncol(test_data)], tau = opt.tau.2, num_iter = opt.k.e2,
                        num_balls = 2)
          result <- as.factor(model2$final_predictions)
          ecccds.auc[2] <- auc(as.numeric(as.factor(result)), as.numeric(as.factor(test_data$Class)))

          # ACDC k = 3
          model3 <- acdc(x = train_data[,-ncol(train_data)], y = train_data$Class,
                        test_data = test_data[,-ncol(test_data)], tau = opt.tau.3, num_iter = opt.k.e3,
                        num_balls = 3)
          result <- as.factor(model3$final_predictions)
          ecccds.auc[3] <- auc(as.numeric(as.factor(result)), as.numeric(as.factor(test_data$Class)))

          # ACDC k = 4
          model4 <- acdc(x = train_data[,-ncol(train_data)], y = train_data$Class,
                        test_data = test_data[,-ncol(test_data)], tau = opt.tau.4, num_iter = opt.k.e4,
                        num_balls = 4)
          result <- as.factor(model4$final_predictions)
          ecccds.auc[4] <- auc(as.numeric(as.factor(result)), as.numeric(as.factor(test_data$Class)))

          # ACDC k = 5
          model5 <- acdc(x = train_data[,-ncol(train_data)], y = train_data$Class,
                        test_data = test_data[,-ncol(test_data)], tau = opt.tau.5, num_iter = opt.k.e5,
                        num_balls = 5)
          result <- as.factor(model5$final_predictions)
          ecccds.auc[5] <- auc(as.numeric(as.factor(result)), as.numeric(as.factor(test_data$Class)))

          max.auc <- c(cccd.auc,ecccds.auc)
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
}
