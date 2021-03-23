gcn_in_out <- function(max_atoms, feature_dim,
    gcn_units, gcn_activation,
    fc_units, fc_activation) {
    if (length(unique(c(length(gcn_units), length(gcn_activation)))) != 1) {
        stop("the number of layers for gcn should be the same")
    }
    if (length(unique(c(length(fc_units), length(fc_activation)))) != 1) {
        stop("the number of layers for fc in gcn_in_out should be the same")
    }
    
    if (!exists("Layer")) {
        # from keras library
        Layer <- function(classname, initialize, build = NULL, call = NULL, 
            compute_output_shape = NULL, ...,  
            inherit = tensorflow::tf$keras$layers$Layer) {
            defs <- list(
                initialize = initialize,
                build = build,
                call = call,
                compute_output_shape = compute_output_shape
            )
            defs <- Filter(Negate(is.null), defs)
            defs <- append(defs, list(...))
            
            
            # allow using the initialize method
            if ("initialize" %in% names(defs)) {
                if (!is.null(defs$`__init__`))
                    stop("You should not specify both __init__ and initialize methods.", call.=FALSE)
                
                defs[["__init__"]] <- defs$initialize
            }
            
            # automatically add the `self` argument
            defs <- lapply(defs, function(x) {
                
                if (inherits(x, "function")) {
                    formals(x) <- append(
                        pairlist(self = NULL),
                        formals(x)
                    )
                }
                
                x
            })
            
            # makes the function return NULL. `__init__` in python must always return None
            defs$`__init__` <- wrap_return_null(defs$`__init__`)
            
            # allow inheriting from custom created layers
            if (!is.null(attr(inherit, "layer")))
                inherit <- attr(inherit, "layer")
            
            layer <- reticulate::PyClass(
                classname = classname,
                defs = defs,
                inherit = inherit
            )
            
            # build the function to be used
            f <- function() {
                .args <- as.list(match.call())[-c(1)]
                .args <- .args[names(.args) != "object"]
                create_layer(layer, object, .args)
            }
            formals(f) <- append(
                list(object = quote(expr=)),
                formals(initialize)
            )
            attr(f, "layer") <- layer
            f
        }
        
        # makes the function return NULL. `__init__` in python must always return None.
        wrap_return_null <- function(f) {
            body(f)[[length(body(f)) + 1]] <- substitute(return(NULL))
            f
        }
    }
    
    layer_multi_linear <- Layer(
        classname <- "MultiLinear", 
        initialize <- function(units, ...) {
            super()$`__init__`(...)
            self$units <- units
        },
        build <- function(input_shape) {
            self$kernel <- self$add_weight(shape = shape(input_shape[3], self$units))
            self$bias <- self$add_weight(shape = shape(self$units))
        },
        call <- function(inputs, ...) {
            keras::k_dot(inputs, self$kernel) + self$bias
        },
        get_config <- function() {
            list(
                units <- self$units,
                name <- self$name
            )
        }
    )
    
    inputA <- keras::layer_input(shape = c(max_atoms, max_atoms))
    inputX <- keras::layer_input(shape = c(max_atoms, feature_dim))
    x <- inputX
    
    for (i in seq_len(length(gcn_units))) {
        assign("temp_units", gcn_units[i], envir = globalenv())
        x <- keras::layer_dot(c(inputA, x), axes = 0)
        x  <- x %>% 
            layer_multi_linear(units = temp_units) %>%
            keras::layer_activation(activation = gcn_activation[i])
    }
    rm(temp_units, envir = globalenv())
    
    output <- keras::layer_global_average_pooling_1d()(x)
    for (i in seq_len(length(fc_units))) {
        output <- output %>% 
            keras::layer_dense(units = fc_units[i], activation = fc_activation[i])
    }
    
    result <- NULL
    result$inputA <- inputA
    result$inputX <- inputX
    result$output <- output
    result
}



cnn_in_out <- function(length_seq, fingerprint_size,
    embedding_layer = TRUE,
    num_tokens, embedding_dim,
    cnn_filters, cnn_kernel_size, cnn_activation,
    fc_units, fc_activation) {
    if (length(unique(c(length(cnn_filters),
        length(cnn_kernel_size),
        length(cnn_activation)))) != 1) {
        stop("the number of layers for cnn should be the same")
    }
    if (length(unique(c(length(fc_units), length(fc_activation)))) != 1) {
        stop("the number of layers for fc in cnn_in_out should be the same")
    }
    
    if (embedding_layer) {
        input <- keras::layer_input(shape = c(length_seq))
        x <- input %>%
            keras::layer_embedding(input_dim = num_tokens + 1,
                output_dim = embedding_dim,
                input_length = length_seq,
                mask_zero = TRUE)
    } else {
        input <- keras::layer_input(shape = c(fingerprint_size, 1))
        x <- input
    }
    
    for (i in seq_len(length(cnn_filters))) {
        x <- x %>%
            keras::layer_conv_1d(filters = cnn_filters[i],
                kernel_size = cnn_kernel_size[i],
                activation = cnn_activation[i])
    }
    
    output <- x %>% keras::layer_global_max_pooling_1d()
    for (i in seq_len(length(fc_units))) {
        output <- output %>% 
            keras::layer_dense(units = fc_units[i], activation = fc_activation[i])
    }
    
    result <- NULL
    result$input <- input
    result$output <- output
    result
}



rnn_in_out <- function(length_seq, fingerprint_size,
    embedding_layer = TRUE,
    num_tokens, embedding_dim,
    rnn_type, rnn_bidirectional,
    rnn_units, rnn_activation,
    fc_units, fc_activation) {
    if (length(unique(c(length(rnn_type),
        length(rnn_bidirectional),
        length(rnn_units),
        length(rnn_activation)))) != 1) {
        stop("the number of layers for rnn should be the same")
    }
    if (length(unique(c(length(fc_units), length(fc_activation)))) != 1) {
        stop("the number of layers for fc in rnn_in_out should be the same")
    }
    
    if (embedding_layer) {
        input <- keras::layer_input(shape = c(length_seq))
        x <- input %>%
            keras::layer_embedding(input_dim = num_tokens + 1,
                output_dim = embedding_dim,
                input_length = length_seq,
                mask_zero = TRUE)
    } else {
        input <- keras::layer_input(shape = c(fingerprint_size, 1))
        x <- input
    }
    
    rnn_return_sequences <- rep(TRUE, length(rnn_units))
    rnn_return_sequences[length(rnn_units)] <- FALSE
    for (i in seq_len(length(rnn_units))) {
        if (rnn_bidirectional[i] & rnn_type[i] == "lstm") {
            x <- x %>%
                keras::bidirectional(keras::layer_lstm(units = rnn_units[i],
                    activation = rnn_activation[i],
                    return_sequences = rnn_return_sequences[i]))
        }
        if (rnn_bidirectional[i] & rnn_type[i] == "gru") {
            x <- x %>%
                keras::bidirectional(keras::layer_gru(units = rnn_units[i],
                    activation = rnn_activation[i],
                    return_sequences = rnn_return_sequences[i]))
        }
        if (!rnn_bidirectional[i] & rnn_type[i] == "lstm") {
            x <- x %>%
                keras::layer_lstm(units = rnn_units[i],
                    activation = rnn_activation[i],
                    return_sequences = rnn_return_sequences[i])
        }
        if (!rnn_bidirectional[i] & rnn_type[i] == "gru") {
            x <- x %>%
                keras::layer_gru(units = rnn_units[i],
                    activation = rnn_activation[i],
                    return_sequences = rnn_return_sequences[i])
        }
    }
    
    output <- x
    for (i in seq_len(length(fc_units))) {
        output <- output %>% 
            keras::layer_dense(units = fc_units[i], activation = fc_activation[i])
    }
    
    result <- NULL
    result$input <- input
    result$output <- output
    result
}



mlp_in_out <- function(length_seq, fingerprint_size,
    embedding_layer = TRUE,
    num_tokens, embedding_dim,
    fc_units, fc_activation) {
    if (length(unique(c(length(fc_units), length(fc_activation)))) != 1) {
        stop("the number of layers for fc in mlp_in_out should be the same")
    }
    
    if (embedding_layer) {
        input <- keras::layer_input(shape = c(length_seq))
        x <- input %>%
            keras::layer_embedding(input_dim = num_tokens + 1,
                output_dim = embedding_dim,
                input_length = length_seq,
                mask_zero = TRUE)
    } else {
        input <- keras::layer_input(shape = c(fingerprint_size, 1))
        x <- input
    }
    
    output <- keras::layer_flatten(x)
    for (i in seq_len(length(fc_units))) {
        output <- output %>% 
            keras::layer_dense(units = fc_units[i], activation = fc_activation[i])
    }
    
    result <- NULL
    result$input <- input
    result$output <- output
    result
}
