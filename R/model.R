fit_cpi <- function(smiles = NULL, AAseq = NULL, outcome,
    convert_canonical_smiles = TRUE,
    compound_type = NULL, compound_max_atoms,
    compound_length_seq, protein_length_seq,
    compound_embedding_dim, protein_embedding_dim,
    protein_ngram_max = 1, protein_ngram_min = 1,
    smiles_val = NULL, AAseq_val = NULL, outcome_val = NULL,
    net_args = list(
        compound,
        compound_args,
        protein,
        protein_args,
        fc_units = c(1),
        fc_activation = c("linear"), ...),
    net_names = list(
        name_compound_max_atoms = NULL,
        name_compound_feature_dim = NULL,
        name_compound_fingerprint_size = NULL,
        name_compound_embedding_layer = NULL,
        name_compound_length_seq = NULL,
        name_compound_num_tokens = NULL,
        name_compound_embedding_dim = NULL,
        name_protein_length_seq = NULL,
        name_protein_num_tokens = NULL,
        name_protein_embedding_dim = NULL),
    preprocessor_only = FALSE,
    preprocessing = list(
        outcome = NULL,
        outcome_val = NULL,
        convert_canonical_smiles = NULL,
        canonical_smiles = NULL,
        compound_type = NULL,
        compound_max_atoms = NULL,
        compound_A_pad = NULL,
        compound_X_pad = NULL,
        compound_A_pad_val = NULL,
        compound_X_pad_val = NULL,
        compound_fingerprint = NULL,
        compound_fingerprint_val = NULL,
        smiles_encode_pad = NULL,
        smiles_val_encode_pad = NULL,
        compound_lenc = NULL,
        compound_length_seq = NULL,
        compound_num_tokens = NULL,
        compound_embedding_dim = NULL,
        AAseq_encode_pad = NULL,
        AAseq_val_encode_pad = NULL,
        protein_lenc = NULL,
        protein_length_seq = NULL,
        protein_num_tokens = NULL,
        protein_embedding_dim = NULL,
        protein_ngram_max = NULL,
        protein_ngram_min = NULL),
    batch_size, use_generator = FALSE,
    validation_split = 0, ...) {
    name_compound_max_atoms <- ifelse(!is.null(net_names$name_compound_max_atoms),
        net_names$name_compound_max_atoms,
        "max_atoms")
    name_compound_feature_dim <- ifelse(!is.null(net_names$name_compound_feature_dim),
        net_names$name_compound_feature_dim,
        "feature_dim")
    name_compound_fingerprint_size <- ifelse(!is.null(net_names$name_compound_fingerprint_size),
        net_names$name_compound_fingerprint_size,
        "fingerprint_size")
    name_compound_embedding_layer <- ifelse(!is.null(net_names$name_compound_embedding_layer),
        net_names$name_compound_embedding_layer,
        "embedding_layer")
    name_compound_length_seq <- ifelse(!is.null(net_names$name_compound_length_seq),
        net_names$name_compound_length_seq,
        "length_seq")
    name_compound_num_tokens <- ifelse(!is.null(net_names$name_compound_num_tokens),
        net_names$name_compound_num_tokens,
        "num_tokens")
    name_compound_embedding_dim <- ifelse(!is.null(net_names$name_compound_embedding_dim),
        net_names$name_compound_embedding_dim,
        "embedding_dim")
    name_protein_length_seq <- ifelse(!is.null(net_names$name_protein_length_seq),
        net_names$name_protein_length_seq,
        "length_seq")
    name_protein_num_tokens <- ifelse(!is.null(net_names$name_protein_num_tokens),
        net_names$name_protein_num_tokens,
        "num_tokens")
    name_protein_embedding_dim <- ifelse(!is.null(net_names$name_protein_embedding_dim),
        net_names$name_protein_embedding_dim,
        "embedding_dim")
    
    result <- NULL
    A_pad <- NULL
    X_pad <- NULL
    fp <- NULL
    smiles_encode_pad <- NULL
    AAseq_encode_pad <- NULL
    A_pad_val <- NULL
    X_pad_val <- NULL
    fp_val <- NULL
    smiles_val_encode_pad <- NULL
    AAseq_val_encode_pad <- NULL
    result$preprocessing <- NULL
    if (any(unlist(lapply(preprocessing, is.null))) | preprocessor_only) {
        ### check data
        message("checking sequences...")
        checked_seq <- seq_check(smiles, AAseq, outcome)
        smiles <- checked_seq$smiles
        AAseq <- checked_seq$AAseq
        outcome <- checked_seq$outcome
        result$preprocessing$outcome <- outcome
        result$preprocessing$removed_smiles <- checked_seq$removed_smiles
        result$preprocessing$removed_AAseq <- checked_seq$removed_AAseq
        
        if (all(!is.null(smiles_val), !is.null(outcome_val)) |
            all(!is.null(AAseq_val), !is.null(outcome_val))) {
            if (!all(ncol(smiles) == ncol(cbind(smiles_val)),
                ncol(AAseq) == ncol(cbind(AAseq_val)),
                ncol(outcome) == ncol(data.matrix(outcome_val)))) {
                stop("dimension of validation data is not the same to the dimension of data")
            }
            checked_seq_val <- seq_check(smiles_val, AAseq_val, outcome_val)
            smiles_val <- checked_seq_val$smiles
            AAseq_val <- checked_seq_val$AAseq
            outcome_val <- checked_seq_val$outcome
            result$preprocessing$outcome_val <- outcome_val
        }
        
        
        ### compound
        if (is.null(smiles)) {
            net_args$compound = net_args$compound_args <- NULL
        } else {
            message("preprocessing for compounds...")
            preprocessed_compound <- seq_preprocessing(smiles = smiles,
                AAseq = NULL,
                compound_type,
                convert_canonical_smiles,
                compound_max_atoms,
                compound_length_seq,
                lenc = NULL)
            
            if (compound_type == "graph") {
                A_pad <- preprocessed_compound$A_pad
                X_pad <- preprocessed_compound$X_pad
                net_args$compound_args[[name_compound_max_atoms]] <- compound_max_atoms
                net_args$compound_args[[name_compound_feature_dim]] <- dim(X_pad[[1]])[3]
            }
            if (compound_type == "fingerprint") {
                fp <- preprocessed_compound$fp
                net_args$compound_args[[name_compound_fingerprint_size]] <- ncol(fp[[1]])
                net_args$compound_args[[name_compound_embedding_layer]] <- FALSE
            }
            if (compound_type == "sequence") {
                smiles_encode_pad <- preprocessed_compound$sequences_encode_pad
                net_args$compound_args[[name_compound_length_seq]] <- compound_length_seq
                net_args$compound_args[[name_compound_num_tokens]] <- preprocessed_compound$num_tokens
                net_args$compound_args[[name_compound_embedding_dim]] <- compound_embedding_dim
            }
            
            # preprocessing for validation data of compounds
            if (all(!is.null(smiles_val), !is.null(outcome_val))) {
                preprocessed_compound_val <- seq_preprocessing(smiles = smiles_val,
                    AAseq = NULL,
                    compound_type,
                    convert_canonical_smiles,
                    compound_max_atoms,
                    compound_length_seq,
                    lenc = preprocessed_compound$lenc)
                # graph
                result$preprocessing$compound_A_pad_val = A_pad_val <-
                    preprocessed_compound_val$A_pad
                result$preprocessing$compound_X_pad_val = X_pad_val <-
                    preprocessed_compound_val$X_pad
                # fingerprint
                result$preprocessing$compound_fingerprint_val = fp_val <-
                    preprocessed_compound_val$fp
                # sequence
                result$preprocessing$smiles_val_encode_pad = smiles_val_encode_pad <-
                    preprocessed_compound_val$sequences_encode_pad
            }
        }
        
        
        ### protein
        if (is.null(AAseq)) {
            net_args$protein = net_args$protein_args <- NULL
        } else {
            message("preprocessing for proteins...")
            preprocessed_protein <- seq_preprocessing(smiles = NULL,
                AAseq = AAseq,
                type = "sequence",
                convert_canonical_smiles,
                max_atoms,
                protein_length_seq,
                lenc = NULL,
                protein_ngram_max,
                protein_ngram_min)
            
            AAseq_encode_pad <- preprocessed_protein$sequences_encode_pad
            net_args$protein_args[[name_protein_length_seq]] <- protein_length_seq
            net_args$protein_args[[name_protein_num_tokens]] <- preprocessed_protein$num_tokens
            net_args$protein_args[[name_protein_embedding_dim]] <- protein_embedding_dim
            
            # preprocessing for validation data of proteins
            if (all(!is.null(AAseq_val), !is.null(outcome_val))) {
                preprocessed_protein_val <- seq_preprocessing(smiles = NULL,
                    AAseq = AAseq_val,
                    type = "sequence",
                    convert_canonical_smiles,
                    compound_max_atoms,
                    protein_length_seq,
                    lenc = preprocessed_protein$lenc,
                    protein_ngram_max,
                    protein_ngram_min)
                
                result$preprocessing$AAseq_val_encode_pad = AAseq_val_encode_pad <-
                    preprocessed_protein_val$sequences_encode_pad
            }
        }
        
        
        # preprocessing result
        if (!is.null(smiles)) {
            result$preprocessing$canonical_smiles <- preprocessed_compound$canonical_smiles
            result$preprocessing$convert_canonical_smiles <- preprocessed_compound$convert_canonical_smiles
            result$preprocessing$compound_type <- compound_type
            if (compound_type == "graph") {
                result$preprocessing$compound_max_atoms <- compound_max_atoms
                result$preprocessing$compound_A_pad <- A_pad
                result$preprocessing$compound_X_pad <- X_pad
            }
            if (compound_type == "fingerprint") {
                result$preprocessing$compound_fingerprint <- fp
            }
            if (compound_type == "sequence") {
                result$preprocessing$compound_lenc <- preprocessed_compound$lenc
                result$preprocessing$smiles_encode_pad <- smiles_encode_pad
                result$preprocessing$compound_length_seq <- compound_length_seq
                result$preprocessing$compound_num_tokens <-
                    net_args$compound_args[[name_compound_num_tokens]]
                result$preprocessing$compound_embedding_dim <- compound_embedding_dim
            }
        }
        if (!is.null(AAseq)) {
            result$preprocessing$protein_lenc <- preprocessed_protein$lenc
            result$preprocessing$AAseq_encode_pad <- AAseq_encode_pad
            result$preprocessing$protein_length_seq <- protein_length_seq
            result$preprocessing$protein_num_tokens <-
                net_args$protein_args[[name_protein_num_tokens]]
            result$preprocessing$protein_embedding_dim <- protein_embedding_dim
            result$preprocessing$protein_ngram_max <- protein_ngram_max
            result$preprocessing$protein_ngram_min <- protein_ngram_min
        }
        
        if (preprocessor_only) return(result)
        
    } else {
        ### from the preprocessed result
        result$preprocessing <- preprocessing
        outcome <- preprocessing$outcome
        compound_type <- preprocessing$compound_type
        
        # compound
        if (is.null(compound_type)) {
            net_args$compound = net_args$compound_args <- NULL
        } else {
            if (compound_type == "graph") {
                net_args$compound_args[[name_compound_max_atoms]] <- preprocessing$compound_max_atoms
                A_pad <- preprocessing$compound_A_pad
                X_pad <- preprocessing$compound_X_pad
                net_args$compound_args[[name_compound_feature_dim]] <- dim(X_pad[[1]])[3]
            }
            if (compound_type == "fingerprint") {
                fp <- preprocessing$compound_fingerprint
                net_args$compound_args[[name_compound_fingerprint_size]] <- ncol(fp[[1]])
                net_args$compound_args[[name_compound_embedding_layer]] <- FALSE
            }
            if (compound_type == "sequence") {
                net_args$compound_args[[name_compound_length_seq]] <- preprocessing$compound_length_seq
                smiles_encode_pad <- preprocessing$smiles_encode_pad
                net_args$compound_args[[name_compound_num_tokens]] <- preprocessing$compound_num_tokens
                net_args$compound_args[[name_compound_embedding_dim]] <- preprocessing$compound_embedding_dim
            }
        }
        # proteins
        if (is.null(preprocessing$AAseq_encode_pad)) {
            net_args$protein = net_args$protein_args <- NULL
        } else {
            net_args$protein_args[[name_protein_length_seq]] <- preprocessing$protein_length_seq
            AAseq_encode_pad <- preprocessing$AAseq_encode_pad
            net_args$protein_args[[name_protein_num_tokens]] <- preprocessing$protein_num_tokens
            net_args$protein_args[[name_protein_embedding_dim]] <- preprocessing$protein_embedding_dim
        }
        
        # validation data for compounds
        if (all(any(!is.null(preprocessing$compound_A_pad_val),
            !is.null(preprocessing$compound_X_pad_val),
            !is.null(preprocessing$compound_fingerprint_val),
            !is.null(preprocessing$smiles_val_encode_pad)),
            !is.null(preprocessing$outcome_val))) {
            outcome_val <- preprocessing$outcome_val
            if (all(!is.null(preprocessing$compound_A_pad_val),
                !is.null(preprocessing$compound_X_pad_val))) {
                A_pad_val <- preprocessing$compound_A_pad_val
                X_pad_val <- preprocessing$compound_X_pad_val
            }
            if (!is.null(preprocessing$compound_fingerprint_val)) {
                fp_val <- preprocessing$compound_fingerprint_val
            }
            if (!is.null(preprocessing$smiles_val_encode_pad)) {
                smiles_val_encode_pad <- preprocessing$smiles_val_encode_pad
            }
        }
        
        # validation data for proteins
        if (all(!is.null(preprocessing$AAseq_val_encode_pad),
            !is.null(preprocessing$outcome_val))) {
            outcome_val <- preprocessing$outcome_val
            AAseq_val_encode_pad <- preprocessing$AAseq_val_encode_pad
        }
    }
    
    
    message("fitting model...")
    if (length(unique(c(length(net_args$fc_units),
        length(net_args$fc_activation)))) != 1) {
        stop("the number of layers for fc in net_args should be the same")
    }
    
    if (!is.null(compound_type)) {
        if (compound_type == "graph") {
            if (!is.null(net_args$protein)) {
                compound_net <- do.call(net_args$compound, net_args$compound_args)
                protein_net <- do.call(net_args$protein, net_args$protein_args)
                output <- layer_concatenate(c(compound_net$output, protein_net$output))
                for(i in seq_len(length(net_args$fc_units))) {
                    output <- output %>% 
                        keras::layer_dense(units = net_args$fc_units[i], 
                            activation = net_args$fc_activation[i])
                }
                model <- keras::keras_model(list(
                    compound_net$inputA,
                    compound_net$inputX,
                    protein_net$input),
                    output)
            } else {
                if (length(X_pad) == 1) {
                    compound_net <- do.call(net_args$compound, net_args$compound_args)
                    output <- compound_net$output
                    for(i in seq_len(length(net_args$fc_units))) {
                        output <- output %>% 
                            keras::layer_dense(units = net_args$fc_units[i], 
                                activation = net_args$fc_activation[i])
                    }
                    model <- keras::keras_model(list(
                        compound_net$inputA,
                        compound_net$inputX),
                        output)
                } else {
                    compound1_net <- do.call(net_args$compound, net_args$compound_args)
                    compound2_net <- do.call(net_args$compound, net_args$compound_args)
                    output <- layer_concatenate(c(compound1_net$output, compound2_net$output))
                    for(i in seq_len(length(net_args$fc_units))) {
                        output <- output %>% 
                            keras::layer_dense(units = net_args$fc_units[i], 
                                activation = net_args$fc_activation[i])
                    }
                    model <- keras::keras_model(list(
                        compound1_net$inputA,
                        compound1_net$inputX,
                        compound2_net$inputA,
                        compound2_net$inputX),
                        output)
                }
            }
        } else if (compound_type == "fingerprint") {
            if (!is.null(net_args$protein)) {
                compound_net <- do.call(net_args$compound, net_args$compound_args)
                protein_net <- do.call(net_args$protein, net_args$protein_args)
                output <- layer_concatenate(c(compound_net$output, protein_net$output))
                for(i in seq_len(length(net_args$fc_units))) {
                    output <- output %>% 
                        keras::layer_dense(units = net_args$fc_units[i], 
                            activation = net_args$fc_activation[i])
                }
                model <- keras::keras_model(list(compound_net$input, protein_net$input),
                    output)
            } else {
                if (length(fp) == 1) {
                    compound_net <- do.call(net_args$compound, net_args$compound_args)
                    output <- compound_net$output
                    for(i in seq_len(length(net_args$fc_units))) {
                        output <- output %>% 
                            keras::layer_dense(units = net_args$fc_units[i], 
                                activation = net_args$fc_activation[i])
                    }
                    model <- keras::keras_model(compound_net$input, output)
                } else {
                    compound1_net <- do.call(net_args$compound, net_args$compound_args)
                    compound2_net <- do.call(net_args$compound, net_args$compound_args)
                    output <- layer_concatenate(c(compound1_net$output, compound2_net$output))
                    for(i in seq_len(length(net_args$fc_units))) {
                        output <- output %>% 
                            keras::layer_dense(units = net_args$fc_units[i], 
                                activation = net_args$fc_activation[i])
                    }
                    model <- keras::keras_model(list(compound1_net$input, compound2_net$input),
                        output)
                }
            }
        } else if (compound_type == "sequence") {
            if (!is.null(net_args$protein)) {
                compound_net <- do.call(net_args$compound, net_args$compound_args)
                protein_net <- do.call(net_args$protein, net_args$protein_args)
                output <- layer_concatenate(c(compound_net$output, protein_net$output))
                for(i in seq_len(length(net_args$fc_units))) {
                    output <- output %>% 
                        keras::layer_dense(units = net_args$fc_units[i], 
                            activation = net_args$fc_activation[i])
                }
                model <- keras::keras_model(list(compound_net$input, protein_net$input),
                    output)
            } else {
                if (length(smiles_encode_pad) == 1) {
                    compound_net <- do.call(net_args$compound, net_args$compound_args)
                    output <- compound_net$output
                    for(i in seq_len(length(net_args$fc_units))) {
                        output <- output %>% 
                            keras::layer_dense(units = net_args$fc_units[i], 
                                activation = net_args$fc_activation[i])
                    }
                    model <- keras::keras_model(compound_net$input, output)
                } else {
                    compound1_net <- do.call(net_args$compound, net_args$compound_args)
                    compound2_net <- do.call(net_args$compound, net_args$compound_args)
                    output <- layer_concatenate(c(compound1_net$output, compound2_net$output))
                    for(i in seq_len(length(net_args$fc_units))) {
                        output <- output %>% 
                            keras::layer_dense(units = net_args$fc_units[i], 
                                activation = net_args$fc_activation[i])
                    }
                    model <- keras::keras_model(list(compound1_net$input, compound2_net$input),
                        output)
                }
            }
        } else {
            stop("model cannnot be found")
        }
    } else {
        if (length(AAseq_encode_pad) == 1) {
            protein_net <- do.call(net_args$protein, net_args$protein_args)
            output <- protein_net$output
            for(i in seq_len(length(net_args$fc_units))) {
                output <- output %>% 
                    keras::layer_dense(units = net_args$fc_units[i], 
                        activation = net_args$fc_activation[i])
            }
            model <- keras::keras_model(list(protein_net$input), output)
        } else {
            protein1_net <- do.call(net_args$protein, net_args$protein_args)
            protein2_net <- do.call(net_args$protein, net_args$protein_args)
            output <- layer_concatenate(c(protein1_net$output, protein2_net$output))
            for(i in seq_len(length(net_args$fc_units))) {
                output <- output %>% 
                    keras::layer_dense(units = net_args$fc_units[i],
                        activation = net_args$fc_activation[i])
            }
            model <- keras::keras_model(list(protein1_net$input, protein2_net$input), output)
        }
    }
    
    net_args$object <- model
    net_args$compound <- NULL
    net_args$compound_args <- NULL
    net_args$protein <- NULL
    net_args$protein_args <- NULL
    net_args$fc_units <- NULL
    net_args$fc_activation <- NULL
    model <- do.call(keras::compile, net_args)
    
    
    validation_data <- NULL
    validation_steps <- NULL
    if (!use_generator) {
        ### without generator
        if (all(any(!is.null(A_pad_val),
            !is.null(X_pad_val),
            !is.null(fp_val),
            !is.null(smiles_val_encode_pad),
            !is.null(AAseq_val_encode_pad)),
            !is.null(outcome_val))) {
            if (exists("fp_val")) {
                for (i in seq_len(length(fp_val))) {
                    dim(fp_val[[i]]) <- c(nrow(fp_val[[i]]), ncol(fp_val[[i]]), 1)
                }
            }
            for (y in c("A_pad_val[[1]]", "X_pad_val[[1]]", "A_pad_val[[2]]", "X_pad_val[[2]]",
                "fp_val[[1]]", "fp_val[[2]]",
                "smiles_val_encode_pad[[1]]", "smiles_val_encode_pad[[2]]",
                "AAseq_val_encode_pad[[1]]", "AAseq_val_encode_pad[[2]]")) {
                if (!is.null(purrr::pluck(eval(parse(text = strsplit(y, split = "\\[\\[")[[1]][1])),
                    as.integer(gsub("[[:alpha:][:punct:]]", "", y))))) {
                    validation_data <- c(validation_data, list(eval(parse(text = y))))
                }
            }
            validation_data <- list(validation_data, outcome_val)
        }
        x <- NULL
        if (exists("fp")) {
            for (i in seq_len(length(fp))) {
                dim(fp[[i]]) <- c(nrow(fp[[i]]), ncol(fp[[i]]), 1)
            }
        }
        for (y in c("A_pad[[1]]", "X_pad[[1]]", "A_pad[[2]]", "X_pad[[2]]",
            "fp[[1]]", "fp[[2]]",
            "smiles_encode_pad[[1]]", "smiles_encode_pad[[2]]",
            "AAseq_encode_pad[[1]]", "AAseq_encode_pad[[2]]")) {
            if (!is.null(purrr::pluck(eval(parse(text = strsplit(y, split = "\\[\\[")[[1]][1])),
                as.integer(gsub("[[:alpha:][:punct:]]", "", y))))) {
                x <- c(x, list(eval(parse(text = y))))
            }
        }
        
        model %>% fit(x, outcome, batch_size = batch_size,
            validation_split = validation_split,
            validation_data = validation_data, ...)
        sampling_generator <- NULL
    } else {
        ### with generator
        idx <- sample(seq_len(nrow(outcome)))
        train_idx <- seq_len(nrow(outcome)) %in%
            idx[seq_len(round(nrow(outcome) * (1 - validation_split)))]
        
        if (all(any(!is.null(A_pad_val),
            !is.null(X_pad_val),
            !is.null(fp_val),
            !is.null(smiles_val_encode_pad),
            !is.null(AAseq_val_encode_pad)),
            !is.null(outcome_val)) | validation_split) {
            if (all(any(!is.null(A_pad_val),
                !is.null(X_pad_val),
                !is.null(fp_val),
                !is.null(smiles_val_encode_pad),
                !is.null(AAseq_val_encode_pad)),
                !is.null(outcome_val))) {
                ### with validation data
                x <- NULL
                if (exists("fp")) {
                    for (i in seq_len(length(fp))) {
                        dim(fp[[i]]) <- c(nrow(fp[[i]]), ncol(fp[[i]]), 1)
                    }
                }
                for (y in c("A_pad[[1]]", "X_pad[[1]]", "A_pad[[2]]", "X_pad[[2]]",
                    "fp[[1]]", "fp[[2]]",
                    "smiles_encode_pad[[1]]", "smiles_encode_pad[[2]]",
                    "AAseq_encode_pad[[1]]", "AAseq_encode_pad[[2]]")) {
                    if (!is.null(purrr::pluck(eval(parse(text = strsplit(y, split = "\\[\\[")[[1]][1])),
                        as.integer(gsub("[[:alpha:][:punct:]]", "", y))))) {
                        x <- c(x, list(eval(parse(text = y))))
                    }
                }
                train_data <- multiple_sampling_generator(x, outcome, batch_size)
                steps_per_epoch <- ceiling(dim(x[[1]])[1]/batch_size)
                
                # validation
                x_val <- NULL
                if (exists("fp_val")) {
                    for (i in seq_len(length(fp_val))) {
                        dim(fp_val[[i]]) <- c(nrow(fp_val[[i]]), ncol(fp_val[[i]]), 1)
                    }
                }
                for (y in c("A_pad_val[[1]]", "X_pad_val[[1]]", "A_pad_val[[2]]", "X_pad_val[[2]]",
                    "fp_val[[1]]", "fp_val[[2]]",
                    "smiles_val_encode_pad[[1]]", "smiles_val_encode_pad[[2]]",
                    "AAseq_val_encode_pad[[1]]", "AAseq_val_encode_pad[[2]]")) {
                    if (!is.null(purrr::pluck(eval(parse(text = strsplit(y, split = "\\[\\[")[[1]][1])),
                        as.integer(gsub("[[:alpha:][:punct:]]", "", y))))) {
                        x_val <- c(x_val, list(eval(parse(text = y))))
                    }
                }
                validation_data <- multiple_sampling_generator(x_val, outcome_val, batch_size)
                validation_steps <- ceiling(dim(x[[1]])[1]/batch_size)
            } else {
                ### split data
                x <- NULL
                if (exists("fp")) {
                    for (i in seq_len(length(fp))) {
                        dim(fp[[i]]) <- c(nrow(fp[[i]]), ncol(fp[[i]]), 1)
                    }
                }
                for (y in c("A_pad[[1]][train_idx,,]", "X_pad[[1]][train_idx,,]",
                    "A_pad[[2]][train_idx,,]", "X_pad[[2]][train_idx,,]",
                    "fp[[1]][train_idx,,,drop = FALSE]", "fp[[2]][train_idx,,,drop = FALSE]",
                    "smiles_encode_pad[[1]][train_idx,]", "smiles_encode_pad[[2]][train_idx,]",
                    "AAseq_encode_pad[[1]][train_idx,]", "AAseq_encode_pad[[2]][train_idx,]")) {
                    if (!is.null(purrr::pluck(eval(parse(text = strsplit(y, split = "\\[\\[")[[1]][1])),
                        as.integer(gsub("[[:alpha:][:punct:]]", "", y))))) {
                        x <- c(x, list(eval(parse(text = y))))
                    }
                }
                train_data <- multiple_sampling_generator(x, outcome[train_idx,,drop = FALSE],
                    batch_size)
                steps_per_epoch <- ceiling(dim(x[[1]])[1]/batch_size)
                
                # validation
                x_val <- NULL
                for (y in c("A_pad[[1]][!train_idx,,]", "X_pad[[1]][!train_idx,,]",
                    "A_pad[[2]][!train_idx,,]", "X_pad[[2]][!train_idx,,]",
                    "fp[[1]][!train_idx,,,drop = FALSE]", "fp[[2]][!train_idx,,,drop = FALSE]",
                    "smiles_encode_pad[[1]][!train_idx,]", "smiles_encode_pad[[2]][!train_idx,]",
                    "AAseq_encode_pad[[1]][!train_idx,]", "AAseq_encode_pad[[2]][!train_idx,]")) {
                    if (!is.null(purrr::pluck(eval(parse(text = strsplit(y, split = "\\[\\[")[[1]][1])),
                        as.integer(gsub("[[:alpha:][:punct:]]", "", y))))) {
                        x_val <- c(x_val, list(eval(parse(text = y))))
                    }
                }
                validation_data <- multiple_sampling_generator(x_val, outcome[!train_idx,,drop = FALSE],
                    batch_size)
                validation_steps <- ceiling(dim(x[[1]])[1]/batch_size)
            }
        } else {
            ### without validation data
            x <- NULL
            if (exists("fp")) {
                for (i in seq_len(length(fp))) {
                    dim(fp[[i]]) <- c(nrow(fp[[i]]), ncol(fp[[i]]), 1)
                }
            }
            for (y in c("A_pad[[1]]", "X_pad[[1]]", "A_pad[[2]]", "X_pad[[2]]",
                "fp[[1]]", "fp[[2]]",
                "smiles_encode_pad[[1]]", "smiles_encode_pad[[2]]",
                "AAseq_encode_pad[[1]]", "AAseq_encode_pad[[2]]")) {
                if (!is.null(purrr::pluck(eval(parse(text = strsplit(y, split = "\\[\\[")[[1]][1])),
                    as.integer(gsub("[[:alpha:][:punct:]]", "", y))))) {
                    x <- c(x, list(eval(parse(text = y))))
                }
            }
            train_data <- multiple_sampling_generator(x, outcome, batch_size)
            steps_per_epoch <- ceiling(dim(x[[1]])[1]/batch_size)
        }
        
        model %>%
            fit(train_data, validation_data = validation_data,
                steps_per_epoch = steps_per_epoch,
                validation_steps = validation_steps, ...)
    }
    
    result$model <- model
    result$sampling_generator <- multiple_sampling_generator
    result
}



predict_cpi <- function(modelRes, smiles = NULL, AAseq = NULL,
    preprocessing = list(
        canonical_smiles = NULL,
        compound_A_pad = NULL,
        compound_X_pad = NULL,
        compound_fingerprint = NULL,
        smiles_encode_pad = NULL,
        AAseq_encode_pad = NULL),
    use_generator = FALSE,
    batch_size = NULL) {
    result <- NULL
    A_pad <- NULL
    X_pad <- NULL
    fp <- NULL
    smiles_encode_pad <- NULL
    AAseq_encode_pad <- NULL
    result$preprocessing <- NULL
    if (any(unlist(lapply(preprocessing, is.null)))) {
        ### check data
        message("checking sequences...")
        checked_seq <- seq_check(smiles, AAseq)
        smiles <- checked_seq$smiles
        AAseq <- checked_seq$AAseq
        
        
        ### compound
        if (!is.null(smiles)) {
            message("preprocessing for compounds...")
            preprocessed_compound <- seq_preprocessing(smiles = smiles,
                AAseq = NULL,
                modelRes$preprocessing$compound_type,
                modelRes$preprocessing$convert_canonical_smiles,
                modelRes$preprocessing$compound_max_atoms,
                modelRes$preprocessing$compound_length_seq,
                lenc = modelRes$preprocessing$compound_lenc)
            
            if (modelRes$preprocessing$compound_type == "graph") {
                A_pad <- preprocessed_compound$A_pad
                X_pad <- preprocessed_compound$X_pad
                result$preprocessing$compound_A_pad <- A_pad
                result$preprocessing$compound_X_pad <- X_pad
            }
            if (modelRes$preprocessing$compound_type == "fingerprint") {
                fp <- preprocessed_compound$fp
                result$preprocessing$compound_fingerprint <- fp
            }
            if (modelRes$preprocessing$compound_type == "sequence") {
                smiles_encode_pad <- preprocessed_compound$sequences_encode_pad
                result$preprocessing$smiles_encode_pad <- smiles_encode_pad
            }
        }
        
        
        ### protein
        if (!is.null(AAseq)) {
            message("preprocessing for proteins...")
            preprocessed_protein <- seq_preprocessing(smiles = NULL,
                AAseq = AAseq,
                type = "sequence",
                modelRes$preprocessing$convert_canonical_smiles,
                modelRes$preprocessing$max_atoms,
                modelRes$preprocessing$protein_length_seq,
                lenc = modelRes$preprocessing$protein_lenc,
                modelRes$preprocessing$protein_ngram_max,
                modelRes$preprocessing$protein_ngram_min)
            AAseq_encode_pad <- preprocessed_protein$sequences_encode_pad
            result$preprocessing$AAseq_encode_pad <- AAseq_encode_pad
        }
    } else {
        ### from the preprocessed result
        result$preprocessing <- preprocessing
        if (modelRes$preprocessing$compound_type == "graph") {
            A_pad <- preprocessing$compound_A_pad
            X_pad <- preprocessing$compound_X_pad
        }
        if (modelRes$preprocessing$compound_type == "fingerprint") {
            fp <- preprocessing$compound_fingerprint
        }
        if (modelRes$preprocessing$compound_type == "sequence") {
            smiles_encode_pad <- preprocessing$smiles_encode_pad
        }
        AAseq_encode_pad <- preprocessing$AAseq_encode_pad
    }
    
    
    message("predicting model...")
    x <- NULL
    if (exists("fp")) {
        for (i in seq_len(length(fp))) {
            dim(fp[[i]]) <- c(nrow(fp[[i]]), ncol(fp[[i]]), 1)
        }
    }
    for (y in c("A_pad[[1]]", "X_pad[[1]]", "A_pad[[2]]", "X_pad[[2]]",
        "fp[[1]]", "fp[[2]]",
        "smiles_encode_pad[[1]]", "smiles_encode_pad[[2]]",
        "AAseq_encode_pad[[1]]", "AAseq_encode_pad[[2]]")) {
        if (!is.null(purrr::pluck(eval(parse(text = strsplit(y, split = "\\[\\[")[[1]][1])),
            as.integer(gsub("[[:alpha:][:punct:]]", "", y))))) {
            x <- c(x, list(eval(parse(text = y))))
        }
    }
    
    if (!use_generator) {
        values <- predict(modelRes$model, x)
    } else {
        steps <- ceiling(dim(x[[1]])[1]/batch_size)
        values <- predict(modelRes$model,
            modelRes$sampling_generator(
            x, batch_size = batch_size,
            shuffle = FALSE),
            steps = steps)
    }
    result$values <- values
    result
}
