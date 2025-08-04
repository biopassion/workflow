def trim_text(text, max_length):
    if len(text) > max_length:
        return text[:max_length].rsplit(' ', 1)[0] + "..."  # Trim to max_length, avoid cutting words in half
    else:
        return text