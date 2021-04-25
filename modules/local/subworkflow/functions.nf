/*
 * -----------------------------------------------------
 *  Utility functions used in nf-core DSL2 module files
 * -----------------------------------------------------
 */

/*
 * Function to initialise default values and to generate a Groovy Map of available options for nf-core modules
 */
def initOptions(Map args) {
    def Map options = [:]
    args.each{ key, val ->
      options[key] = val ?: [:]
    }
    return options
}
