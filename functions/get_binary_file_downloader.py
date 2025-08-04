# code from 
# https://discuss.streamlit.io/t/how-to-download-file-in-streamlit/1806/27
import os
import base64
def get_binary_file_downloader(bin_file, file_label='File'):
    with open(bin_file, 'rb') as f:
        data = f.read()
    bin_str = base64.b64encode(data).decode()
    href = f'<a href="data:application/octet-stream;base64,{bin_str}" download="{os.path.basename(bin_file)}">Download {file_label}</a>'
    return href