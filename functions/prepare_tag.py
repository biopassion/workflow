import datetime

# generate a tag string 
def prepare_tag(key_word_no_space):
    # time 
    now=datetime.datetime.now()
    now_str= f"{now:%Y_%m_%d}"
    
    return f"{key_word_no_space}_{now_str}"