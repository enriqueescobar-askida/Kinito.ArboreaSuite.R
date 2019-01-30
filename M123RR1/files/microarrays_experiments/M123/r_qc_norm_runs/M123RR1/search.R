# Sébastien Boisvert <s@member.fsf.org>
#
# Array.includes array, object
#
#





Array.recursive_binary_search <- function(array_o, object) {
        first_index = 1
        last_index = length(array_o)

        if(length(array_o)==1){               # the length() == 1 case
                if(array_o[1]==object){
                    return (TRUE)
                }else{
                        return(FALSE)
                }
        }

        middle_index = floor((last_index) / 2)
        middle_value = array_o[middle_index]

        if(object == middle_value) {      # got it!
                return ( TRUE)
        }else if(object > middle_value) {     # search after the middle
                return  (Array.recursive_binary_search(array_o[(middle_index+1):last_index], object))
        }else if(object < middle_value) {     # search before the middle
                return (Array.recursive_binary_search(array_o[1:(middle_index-1)], object))
        }
}

#########################################################

Array.includes <- local(function(array_object, element) {

        array_length <- length(array_object)
        if ( array_length == 0 ) {
                return (FALSE)
        }

        sorted_array <- sort(array_object)
       
        return (Array.recursive_binary_search(sorted_array, element))
})
