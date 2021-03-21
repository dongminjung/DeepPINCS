metric_concordance_index <- function(y_true, y_pred) {
    g <- tensorflow::tf$subtract(
        tensorflow::tf$expand_dims(y_pred, as.integer(-1)), y_pred)
    g <- keras::k_cast(g == 0.0, "float32") * 0.5 + keras::k_cast(g > 0.0, "float32")
    
    f <- tensorflow::tf$subtract(
        tensorflow::tf$expand_dims(y_true, as.integer(-1)), y_true) > 0.0
    f <- tensorflow::tf$linalg$band_part(
        keras::k_cast(f, "float32"), as.integer(-1), as.integer(0))
    
    g <- keras::k_sum(g * f)
    f <- keras::k_sum(f)
    
    tensorflow::tf$where(keras::k_equal(g, 0), 0.0, g/f)
}



metric_f1_score <- function(y_true, y_pred) {
    y_pred <- keras::k_round(y_pred)
    TP <- keras::k_sum(y_pred * y_true)
    precision <- tensorflow::tf$where(keras::k_equal(TP, 0), 0.0, TP / keras::k_sum(y_pred))
    recall <- tensorflow::tf$where(keras::k_equal(TP, 0), 0.0, TP / keras::k_sum(y_true))
    tensorflow::tf$where(keras::k_equal(2 * precision * recall, 0),
        0.0, (2 * precision * recall) / (precision + recall))
}



get_canonical_smiles <- function(smiles) {
    m <- rcdk::parse.smiles(smiles)
    canonical_smiles <- lapply(m, function(x) 
        rcdk::get.smiles(x, flavor = rcdk::smiles.flavors(c("Canonical"))))
    as.vector(unlist(canonical_smiles))
}



get_fingerprint <- function(smiles, ...) {
    m <- rcdk::parse.smiles(smiles)
    fpm <- lapply(m, function(x) rcdk::get.fingerprint(x, ...))
    for (i in seq_len(length(fpm))) {
        if (i == 1) fp <- matrix(0, length(m), fpm[[i]]@nbit)
        fp[i, fpm[[i]]@bits] <- 1
    }
    fp
}



get_graph_structure_node_feature <- function(
    smiles, max_atoms,
    element_list = c(
        "C", "N", "O", "S", "F", "Si", "P", "Cl",
        "Br", "Mg", "Na", "Ca", "Fe", "Al", "I",
        "B", "K", "Se", "Zn", "H", "Cu", "Mn")) {
    A <- list()
    X <- list()
    m <- rcdk::parse.smiles(smiles)
    
    for (i in seq_len(length(m))) {
        adj <- rcdk::get.adjacency.matrix(m[[i]])
        a <- adj + diag(1, dim(adj))
        degree <- rowSums(adj)
        A[[i]] <- diag((degree+1)^(-1/2)) %*% a %*% diag((degree+1)^(-1/2))
        
        atoms <- rcdk::get.atoms(m[[i]])
        num_atoms <- length(degree)
        atomic_symbol <- matrix(0, num_atoms, length(element_list))
        hydrogen_count <- NULL
        for (j in seq_len(num_atoms)) {
            atomic_symbol[j, match(rcdk::get.symbol(atoms[[j]]), element_list)] <- 1
            hydrogen_count <- c(hydrogen_count, rcdk::get.hydrogen.count(atoms[[j]]))
        }
        X[[i]] <- cbind(atomic_symbol, hydrogen_count, degree)
    }
    
    A_pad <- array(0, dim = c(length(A), max_atoms, max_atoms))
    X_pad <- array(0, dim = c(length(X), max_atoms, ncol(purrr::chuck(X, 1))))
    for (i in seq_len(length(m))) {
        n <- dim(A[[i]])[1]
        if (n < max_atoms) {
            A_pad[i,,] <- matlab::padarray(A[[i]], c(max_atoms-n, max_atoms-n), 0, "post")
            X_pad[i,,] <- matlab::padarray(X[[i]], max_atoms-n, 0, "post")
        } else {
            A_pad[i,,] <- A[[i]][seq(max_atoms), seq(max_atoms)]
            X_pad[i,,] <- X[[i]][seq(max_atoms), ]
        }
    }
    
    result <- NULL
    result$A_pad <- A_pad
    result$X_pad <- X_pad
    result$element_list <- element_list
    result$feature_dim <- dim(X_pad)[3]
    result
}



get_seq_encode_pad <- function(sequences, length_seq,
    ngram_max = 1, ngram_min = 1,
    lenc = NULL) {
    subsequences <- NULL
    for (i in ngram_min:ngram_max) {
        subsequences <- cbind(subsequences, stringdist::qgrams(sequences, q = i))
    }
    subsequences <- subsequences[, order(subsequences, decreasing = TRUE), drop = FALSE]
    
    sequences_token <- tokenizers::tokenize_ngrams(gsub("", " ", sequences),
        n = ngram_max,
        n_min = ngram_min,
        lowercase = FALSE,
        ngram_delim = "")
    if (is.null(lenc)) {
        lenc <- CatEncoders::LabelEncoder.fit(colnames(subsequences))
    }
    sequences_encode <- lapply(sequences_token,
        function(x) na.omit(CatEncoders::transform(lenc, x)))
    
    sequences_encode_pad <- sequences_encode %>% 
        keras::pad_sequences(maxlen = length_seq)
    
    result <- NULL
    result$sequences_encode_pad <- sequences_encode_pad
    result$lenc <- lenc
    result$num_tokens <- length(subsequences)
    result
}



seq_check <- function(smiles = NULL, AAseq = NULL, outcome = NULL) {
    result <- NULL
    outcome <- if (!is.null(outcome)) data.matrix(outcome)
    smiles <- if (!is.null(smiles)) cbind(smiles)
    AAseq <- if (!is.null(AAseq)) cbind(AAseq)
    if ((ifelse(is.null(smiles), 0, ncol(smiles)) + 
        ifelse(is.null(AAseq), 0, ncol(AAseq))) > 2) {
        stop("data is not single or paired sequences")
    }
    if (!is.null(outcome)) {
        if (length(unique(c(ifelse(is.null(smiles), nrow(outcome), nrow(smiles)),
            ifelse(is.null(AAseq), nrow(outcome), nrow(AAseq)),
            nrow(outcome)))) != 1) {
            stop("check the number of samples")
        }
    }
    if (!is.null(smiles)) {
        smiles_check <- NULL
        for (i in seq_len(ncol(smiles))) {
            smiles_check[[i]] <- vapply(smiles[,i], function(x) webchem::is.smiles(x, verbose = FALSE),
                TRUE, USE.NAMES = FALSE)
        }
    }
    if (!is.null(AAseq)) {
        AAseq_check <- NULL
        for (i in seq_len(ncol(AAseq))) {
            AAseq_check[[i]] <- grepl("^[A-Z]+$", AAseq[,i])
        }
    }
    if (!is.null(smiles) & !is.null(AAseq)) {
        smiles_AAseq_check <- smiles_check[[1]] & AAseq_check[[1]]
    } else if (!is.null(smiles)) {
        smiles_AAseq_check <- if (ncol(smiles) == 1) smiles_check[[1]] else smiles_check[[1]] & smiles_check[[2]]
    } else {
        smiles_AAseq_check <- if (ncol(AAseq) == 1) AAseq_check[[1]] else AAseq_check[[1]] & AAseq_check[[2]]
    }
    
    if (!all(smiles_AAseq_check)) {
        smiles <- if (!is.null(smiles)) cbind(smiles[smiles_AAseq_check,])
        AAseq <- if (!is.null(AAseq)) cbind(AAseq[smiles_AAseq_check,])
        outcome <- if (!is.null(outcome)) outcome[smiles_AAseq_check,]
        if (!is.null(smiles)) {
            for (i in seq_len(ncol(smiles))) {
                if (!all(smiles_check[[i]])) {
                    result$removed_smiles[[i]] <- which(!smiles_check[[i]])
                    message("at least one of compound sequences may not be valid")
                }
            }
        }
        if (!is.null(AAseq)) {
            for (i in seq_len(ncol(AAseq))) {
                if (!all(AAseq_check[[i]])) {
                    result$removed_AAseq[[i]] <- which(!AAseq_check[[i]])
                    message("at least one of protein sequences may not be valid")
                }
            }
        }
    }
    result$smiles <- smiles
    result$AAseq <- AAseq
    result$outcome <- outcome
    result
}



seq_preprocessing <- function(smiles = NULL,
    AAseq = NULL,
    type,
    convert_canonical_smiles,
    max_atoms,
    length_seq,
    lenc = NULL,
    ngram_max = 1,
    ngram_min = 1) {
    result <- NULL
    if (!is.null(smiles)) {
        if (convert_canonical_smiles) {
            for (i in seq_len(ncol(smiles))) {
                smiles[,i] <- get_canonical_smiles(smiles[,i])
            }
            result$canonical_smiles <- smiles
            result$convert_canonical_smiles <- convert_canonical_smiles
        } else {
            result$convert_canonical_smiles <- convert_canonical_smiles
        }
        sequences <- smiles
    } else {
        sequences <- AAseq
    }
    
    if (type == "graph") {
        A_pad <- NULL
        X_pad <- NULL
        for (i in seq_len(ncol(sequences))) {
            graph_structure_node_feature <-
                get_graph_structure_node_feature(sequences[,i], max_atoms)
            A_pad[[i]] <- graph_structure_node_feature$A_pad
            X_pad[[i]] <- graph_structure_node_feature$X_pad
        }
        result$A_pad <- A_pad
        result$X_pad <- X_pad
    }
    if (type == "fingerprint") {
        fp <- NULL
        for (i in seq_len(ncol(sequences))) {
            fp[[i]] <- get_fingerprint(sequences[,i])
        }
        result$fp <- fp
    }
    if (type == "sequence") {
        if (is.null(lenc)) {
            seq_temp <- sequences
            if (ncol(sequences) == 2) seq_temp <- c(sequences[,1], sequences[,2])
            seq_encode_pad <- get_seq_encode_pad(seq_temp, length_seq,
                ngram_max, ngram_min)
            result$lenc <- seq_encode_pad$lenc
            result$num_tokens <- seq_encode_pad$num_tokens
        }
        sequences_encode_pad <- NULL
        for (i in seq_len(ncol(sequences))) {
            seq_encode_pad <- get_seq_encode_pad(sequences[,i], length_seq,
                ngram_max, ngram_min,
                lenc = lenc)
            sequences_encode_pad[[i]] <- seq_encode_pad$sequences_encode_pad
        }
        result$sequences_encode_pad <- sequences_encode_pad
    }
    result
}



multiple_sampling_generator <- function(X_data, Y_data = NULL, batch_size,
    shuffle = TRUE) {
    if (!shuffle) {
        batch_n <- dim(X_data[[1]])[1]
        assign("batch_start", 0, envir = globalenv())
    }
    
    function() {
        if (!shuffle) {
            rows_to_read <- ifelse(batch_start + batch_size > batch_n,
                batch_n - batch_start, batch_size)
            rows <- (batch_start + 1):(batch_start + rows_to_read)
            new_start <- ifelse(batch_start + batch_size < batch_n,
                batch_start + batch_size, 0)
            if (new_start == 0) {
                if (exists("batch_start")) rm(batch_start, envir = globalenv())
            }
            assign("batch_start", new_start, envir = globalenv())
        } else {
            rows <- sample(seq_len(dim(X_data[[1]])[1]), batch_size, replace = TRUE)
        }
        
        gen_X_data <- NULL
        for (i in seq_len(length(X_data))) {
            temp_X <- eval(parse(text = paste("X_data[[",  i, "]][rows",
                paste(rep(",", length(dim(X_data[[i]]))), collapse = ""),
                "drop = FALSE]", sep = "")))
            gen_X_data <- c(gen_X_data, list(temp_X))
        }
        if (is.null(Y_data)) {
            list(gen_X_data)
        } else {
            gen_Y_data <- eval(parse(text = paste("Y_data[rows",
                paste(rep(",", length(dim(Y_data))), collapse = ""),
                "drop = FALSE]", sep = "")))
            list(gen_X_data, gen_Y_data)
        }
    }
}

