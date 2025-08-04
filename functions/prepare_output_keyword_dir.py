import os

# generate output keyword dir
def prepare_output_keyword_dir(main_dir, key_word_no_space, output_string="output"):
    # output dir 
    output_main_dir = os.path.join(main_dir,output_string)

    # save output dir 
    output_keyword_dir= os.path.join(output_main_dir, key_word_no_space)

    # create dirs if not existed 
    list_dirs = [output_main_dir,output_keyword_dir]
    for directory in list_dirs:
        if not os.path.exists(directory):
            os.makedirs(directory)
    return output_keyword_dir

