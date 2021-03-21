compound_max_atoms <- 50
protein_embedding_dim <- 16
protein_length_seq <- 100
gcn_cnn_cpi <- fit_cpi(smiles = example_cpi[1:100, 1],
                       AAseq = example_cpi[1:100, 2],
                       outcome = example_cpi[1:100, 3],
                       compound_type = "graph",
                       compound_max_atoms = compound_max_atoms,
                       protein_length_seq = protein_length_seq,
                       protein_embedding_dim = protein_embedding_dim,
                       smiles_val = example_cpi[101:200, 1],
                       AAseq_val = example_cpi[101:200, 2],
                       outcome_val = example_cpi[101:200, 3],
                       net_args =
                         list(compound = "gcn_in_out",
                              compound_args = list(gcn_units = c(128, 64),
                                                   gcn_activation = c("relu", "relu"),
                                                   fc_units = c(10),
                                                   fc_activation = c("relu")),
                              protein = "cnn_in_out",
                              protein_args = list(cnn_filters = c(32),
                                                  cnn_kernel_size = c(3),
                                                  cnn_activation = c("relu"),
                                                  fc_units = c(10),
                                                  fc_activation = c("relu")),
                              fc_units = c(1),
                              fc_activation = c("sigmoid"),
                              loss = "binary_crossentropy",
                              optimizer = keras::optimizer_adam(),
                              metrics = "accuracy"),
                       epochs = 1, batch_size = 16,
                       use_generator = TRUE)



test_that("fit_cpi: match the names of arguments", {
  gcn_in_out2 <<- function(max_Atoms, feature_dim, gcn_units, gcn_activation, 
                           fc_units, fc_activation) 
  {
    if (length(unique(c(length(gcn_units), length(gcn_activation)))) != 
        1) {
      stop("the number of layers for gcn should be the same")
    }
    if (length(unique(c(length(fc_units), length(fc_activation)))) != 
        1) {
      stop("the number of layers for fc in gcn_in_out should be the same")
    }
    layer_multi_linear <- keras::Layer(classname <- "MultiLinear", 
                                       initialize <- function(units, ...) {
                                         super()$`__init__`(...)
                                         self$units <- units
                                       }, build <- function(input_shape) {
                                         self$kernel <- self$add_weight(shape = shape(input_shape[3], 
                                                                                      self$units))
                                         self$bias <- self$add_weight(shape = shape(self$units))
                                       }, call <- function(inputs, ...) {
                                         keras::k_dot(inputs, self$kernel) + self$bias
                                       }, get_config <- function() {
                                         list(units <- self$units, name <- self$name)
                                       })
    inputA <- keras::layer_input(shape = c(max_Atoms, max_Atoms))
    inputX <- keras::layer_input(shape = c(max_Atoms, feature_dim))
    x <- inputX
    for (i in seq_len(length(gcn_units))) {
      assign("temp_units", gcn_units[i], envir = globalenv())
      x <- keras::layer_dot(c(inputA, x), axes = 0)
      x <- x %>% layer_multi_linear(units = temp_units) %>% 
        keras::layer_activation(activation = gcn_activation[i])
    }
    rm(temp_units, envir = globalenv())
    output <- (keras::layer_global_average_pooling_1d())(x)
    for (i in seq_len(length(fc_units))) {
      output <- output %>% keras::layer_dense(units = fc_units[i], 
                                              activation = fc_activation[i])
    }
    result <- NULL
    result$inputA <- inputA
    result$inputX <- inputX
    result$output <- output
    result
  }
  
  compound_max_atoms <- 50
  protein_embedding_dim <- 16
  protein_length_seq <- 100
  gcn2_cnn_cpi <- fit_cpi(smiles = example_cpi[1:100, 1],
                          AAseq = example_cpi[1:100, 2],
                          outcome = example_cpi[1:100, 3],
                          compound_type = "graph",
                          compound_max_atoms = compound_max_atoms,
                          protein_length_seq = protein_length_seq,
                          protein_embedding_dim = protein_embedding_dim,
                          smiles_val = example_cpi[101:200, 1],
                          AAseq_val = example_cpi[101:200, 2],
                          outcome_val = example_cpi[101:200, 3],
                          net_args =
                            list(compound = "gcn_in_out2",
                                 compound_args = list(gcn_units = c(128, 64),
                                                      gcn_activation = c("relu", "relu"),
                                                      fc_units = c(10),
                                                      fc_activation = c("relu")),
                                 protein = "cnn_in_out",
                                 protein_args = list(cnn_filters = c(32),
                                                     cnn_kernel_size = c(3),
                                                     cnn_activation = c("relu"),
                                                     fc_units = c(10),
                                                     fc_activation = c("relu")),
                                 fc_units = c(1),
                                 fc_activation = c("sigmoid"),
                                 loss = "binary_crossentropy",
                                 optimizer = keras::optimizer_adam(),
                                 metrics = "accuracy"),
                          net_names = list(name_compound_max_atoms = "max_Atoms"),
                          epochs = 1, batch_size = 16,
                          use_generator = TRUE)
  expect_type(gcn2_cnn_cpi, "list")
})



test_that("fit_cpi: add dropout layer", {
  gcn_in_out3 <<- function(max_atoms, feature_dim, gcn_units, gcn_activation, 
                           fc_units, fc_activation) 
  {
    if (length(unique(c(length(gcn_units), length(gcn_activation)))) != 
        1) {
      stop("the number of layers for gcn should be the same")
    }
    if (length(unique(c(length(fc_units), length(fc_activation)))) != 
        1) {
      stop("the number of layers for fc in gcn_in_out should be the same")
    }
    layer_multi_linear <- keras::Layer(classname <- "MultiLinear", 
                                       initialize <- function(units, ...) {
                                         super()$`__init__`(...)
                                         self$units <- units
                                       }, build <- function(input_shape) {
                                         self$kernel <- self$add_weight(shape = shape(input_shape[3], 
                                                                                      self$units))
                                         self$bias <- self$add_weight(shape = shape(self$units))
                                       }, call <- function(inputs, ...) {
                                         keras::k_dot(inputs, self$kernel) + self$bias
                                       }, get_config <- function() {
                                         list(units <- self$units, name <- self$name)
                                       })
    inputA <- keras::layer_input(shape = c(max_atoms, max_atoms))
    inputX <- keras::layer_input(shape = c(max_atoms, feature_dim))
    x <- inputX
    for (i in seq_len(length(gcn_units))) {
      assign("temp_units", gcn_units[i], envir = globalenv())
      x <- keras::layer_dot(c(inputA, x), axes = 0)
      x <- x %>% layer_multi_linear(units = temp_units) %>% 
        keras::layer_activation(activation = gcn_activation[i]) %>%
        keras::layer_dropout(0.5)
    }
    rm(temp_units, envir = globalenv())
    output <- (keras::layer_global_average_pooling_1d())(x)
    for (i in seq_len(length(fc_units))) {
      output <- output %>% keras::layer_dense(units = fc_units[i], 
                                              activation = fc_activation[i])
    }
    result <- NULL
    result$inputA <- inputA
    result$inputX <- inputX
    result$output <- output
    result
  }
  
  compound_max_atoms <- 50
  protein_embedding_dim <- 16
  protein_length_seq <- 100
  gcn3_cnn_cpi <- fit_cpi(smiles = example_cpi[1:100, 1],
                          AAseq = example_cpi[1:100, 2],
                          outcome = example_cpi[1:100, 3],
                          compound_type = "graph",
                          compound_max_atoms = compound_max_atoms,
                          protein_length_seq = protein_length_seq,
                          protein_embedding_dim = protein_embedding_dim,
                          smiles_val = example_cpi[101:200, 1],
                          AAseq_val = example_cpi[101:200, 2],
                          outcome_val = example_cpi[101:200, 3],
                          net_args =
                            list(compound = "gcn_in_out3",
                                 compound_args = list(gcn_units = c(128, 64),
                                                      gcn_activation = c("relu", "relu"),
                                                      fc_units = c(10),
                                                      fc_activation = c("relu")),
                                 protein = "cnn_in_out",
                                 protein_args = list(cnn_filters = c(32),
                                                     cnn_kernel_size = c(3),
                                                     cnn_activation = c("relu"),
                                                     fc_units = c(10),
                                                     fc_activation = c("relu")),
                                 fc_units = c(1),
                                 fc_activation = c("sigmoid"),
                                 loss = "binary_crossentropy",
                                 optimizer = keras::optimizer_adam(),
                                 metrics = "accuracy"),
                          epochs = 1, batch_size = 16,
                          use_generator = TRUE)
  expect_type(gcn3_cnn_cpi, "list")
})



test_that("fit_cpi: miss gcn_units and gcn_activation", {
  compound_max_atoms <- 50
  protein_embedding_dim <- 16
  protein_length_seq <- 100
  expect_error(
    fit_cpi(smiles = example_cpi[1:100, 1],
            AAseq = example_cpi[1:100, 2],
            outcome = example_cpi[1:100, 3],
            compound_type = "graph",
            compound_max_atoms = compound_max_atoms,
            protein_length_seq = protein_length_seq,
            protein_embedding_dim = protein_embedding_dim,
            smiles_val = example_cpi[101:200, 1],
            AAseq_val = example_cpi[101:200, 2],
            outcome_val = example_cpi[101:200, 3],
            net_args =
              list(compound = "gcn_in_out",
                   compound_args = list(fc_units = c(10),
                                        fc_activation = c("relu")),
                   protein = "cnn_in_out",
                   protein_args = list(cnn_filters = c(32),
                                       cnn_kernel_size = c(3),
                                       cnn_activation = c("relu"),
                                       fc_units = c(10),
                                       fc_activation = c("relu")),
                   fc_units = c(1),
                   fc_activation = c("sigmoid"),
                   loss = "binary_crossentropy",
                   optimizer = keras::optimizer_adam(),
                   metrics = "accuracy"),
            epochs = 1, batch_size = 16,
            use_generator = TRUE))
})



test_that("fit_cpi: incorrect cnn_filters", {
  compound_max_atoms <- 50
  protein_embedding_dim <- 16
  protein_length_seq <- 100
  expect_error(
    fit_cpi(smiles = example_cpi[1:100, 1],
            AAseq = example_cpi[1:100, 2],
            outcome = example_cpi[1:100, 3],
            compound_type = "graph",
            compound_max_atoms = compound_max_atoms,
            protein_length_seq = protein_length_seq,
            protein_embedding_dim = protein_embedding_dim,
            smiles_val = example_cpi[101:200, 1],
            AAseq_val = example_cpi[101:200, 2],
            outcome_val = example_cpi[101:200, 3],
            net_args =
              list(compound = "gcn_in_out",
                   compound_args = list(gcn_units = c(128, 64),
                                        gcn_activation = c("relu", "relu"),
                                        fc_units = c(10),
                                        fc_activation = c("relu")),
                   protein = "cnn_in_out",
                   protein_args = list(cnn_filters = c(32, 32),
                                       cnn_kernel_size = c(3),
                                       cnn_activation = c("relu"),
                                       fc_units = c(10),
                                       fc_activation = c("relu")),
                   fc_units = c(1),
                   fc_activation = c("sigmoid"),
                   loss = "binary_crossentropy",
                   optimizer = keras::optimizer_adam(),
                   metrics = "accuracy"),
            epochs = 1, batch_size = 16,
            use_generator = TRUE))
})



test_that("fit_cpi: incorrect argument of encoder", {
  compound_max_atoms <- 50
  protein_embedding_dim <- 16
  protein_length_seq <- 100
  expect_error(
    fit_cpi(smiles = example_cpi[1:100, 1],
            AAseq = example_cpi[1:100, 2],
            outcome = example_cpi[1:100, 3],
            compound_type = "graph",
            compound_max_atoms = compound_max_atoms,
            protein_length_seq = protein_length_seq,
            protein_embedding_dim = protein_embedding_dim,
            smiles_val = example_cpi[101:200, 1],
            AAseq_val = example_cpi[101:200, 2],
            outcome_val = example_cpi[101:200, 3],
            net_args =
              list(compound = "mlp_in_out",
                   compound_args = list(gcn_units = c(128, 64),
                                        gcn_activation = c("relu", "relu"),
                                        fc_units = c(10),
                                        fc_activation = c("relu")),
                   protein = "cnn_in_out",
                   protein_args = list(cnn_filters = c(32),
                                       cnn_kernel_size = c(3),
                                       cnn_activation = c("relu"),
                                       fc_units = c(10),
                                       fc_activation = c("relu")),
                   fc_units = c(1),
                   fc_activation = c("sigmoid"),
                   loss = "binary_crossentropy",
                   optimizer = keras::optimizer_adam(),
                   metrics = "accuracy"),
            epochs = 1, batch_size = 16,
            use_generator = TRUE))
})



test_that("fit_cpi: incorrect protein encoder", {
  compound_max_atoms <- 50
  protein_embedding_dim <- 16
  protein_length_seq <- 100
  expect_error(
    fit_cpi(smiles = example_cpi[1:100, 1],
            AAseq = example_cpi[1:100, 2],
            outcome = example_cpi[1:100, 3],
            compound_type = "graph",
            compound_max_atoms = compound_max_atoms,
            protein_length_seq = protein_length_seq,
            protein_embedding_dim = protein_embedding_dim,
            smiles_val = example_cpi[101:200, 1],
            AAseq_val = example_cpi[101:200, 2],
            outcome_val = example_cpi[101:200, 3],
            net_args =
              list(compound = "gcn_in_out",
                   compound_args = list(gcn_units = c(128, 64),
                                        gcn_activation = c("relu", "relu"),
                                        fc_units = c(10),
                                        fc_activation = c("relu")),
                   protein = "gcn_in_out",
                   protein_args = list(gcn_units = c(128, 64),
                                       gcn_activation = c("relu", "relu"),
                                       fc_units = c(10),
                                       fc_activation = c("relu")),
                   fc_units = c(1),
                   fc_activation = c("sigmoid"),
                   loss = "binary_crossentropy",
                   optimizer = keras::optimizer_adam(),
                   metrics = "accuracy"),
            epochs = 1, batch_size = 16,
            use_generator = TRUE))
})



test_that("fit_cpi: incorrect input", {
  compound_max_atoms <- 50
  protein_embedding_dim <- 16
  protein_length_seq <- 100
  expect_error(
    fit_cpi(smiles = example_cpi[1:100, 1:2],
            AAseq = example_cpi[1:100, 2],
            outcome = example_cpi[1:100, 3],
            compound_type = "graph",
            compound_max_atoms = compound_max_atoms,
            protein_length_seq = protein_length_seq,
            protein_embedding_dim = protein_embedding_dim,
            smiles_val = example_cpi[101:200, 1:2],
            AAseq_val = example_cpi[101:200, 2],
            outcome_val = example_cpi[101:200, 3],
            net_args =
              list(compound = "gcn_in_out",
                   compound_args = list(gcn_units = c(128, 64),
                                        gcn_activation = c("relu", "relu"),
                                        fc_units = c(10),
                                        fc_activation = c("relu")),
                   protein = "gcn_in_out",
                   protein_args = list(gcn_units = c(128, 64),
                                       gcn_activation = c("relu", "relu"),
                                       fc_units = c(10),
                                       fc_activation = c("relu")),
                   fc_units = c(1),
                   fc_activation = c("sigmoid"),
                   loss = "binary_crossentropy",
                   optimizer = keras::optimizer_adam(),
                   metrics = "accuracy"),
            epochs = 1, batch_size = 16,
            use_generator = TRUE))
})



test_that("fit_cpi: incorrect compound_type", {
  compound_max_atoms <- 50
  protein_embedding_dim <- 16
  protein_length_seq <- 100
  expect_error(
    fit_cpi(smiles = example_cpi[1:100, 1:2],
            AAseq = example_cpi[1:100, 2],
            outcome = example_cpi[1:100, 3],
            compound_type = "sequence",
            compound_max_atoms = compound_max_atoms,
            protein_length_seq = protein_length_seq,
            protein_embedding_dim = protein_embedding_dim,
            smiles_val = example_cpi[101:200, 1:2],
            AAseq_val = example_cpi[101:200, 2],
            outcome_val = example_cpi[101:200, 3],
            net_args =
              list(compound = "gcn_in_out",
                   compound_args = list(gcn_units = c(128, 64),
                                        gcn_activation = c("relu", "relu"),
                                        fc_units = c(10),
                                        fc_activation = c("relu")),
                   protein = "gcn_in_out",
                   protein_args = list(gcn_units = c(128, 64),
                                       gcn_activation = c("relu", "relu"),
                                       fc_units = c(10),
                                       fc_activation = c("relu")),
                   fc_units = c(1),
                   fc_activation = c("sigmoid"),
                   loss = "binary_crossentropy",
                   optimizer = keras::optimizer_adam(),
                   metrics = "accuracy"),
            epochs = 1, batch_size = 16,
            use_generator = TRUE))
})



test_that("predict_cpi: list output", {
  expect_type(
    predict_cpi(gcn_cnn_cpi, example_cpi[201:210, 1], example_cpi[201:210, 2]),
    "list")
})



test_that("predict_cpi: miss batch_size", {
  expect_error(
    predict_cpi(gcn_cnn_cpi, example_cpi[201:210, 1], example_cpi[201:210, 2],
                use_generator = TRUE))
})



test_that("predict_cpi: incorrect input", {
  expect_error(
    predict_cpi(gcn_cnn_cpi, example_cpi[201:210, 1], example_cpi[211:220, 1]))
})



test_that("predict_cpi: incorrect dimension of input", {
  expect_error(
    predict_cpi(gcn_cnn_cpi, example_cpi[201:210, 1], example_cpi[201:209, 1]))
})