def get_config_key(file):
    if "stringtie" in file:
        return "stringtie"
    else:
        return os.path.basename(file)
