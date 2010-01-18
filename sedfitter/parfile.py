def read(filename):
    '''
    Read a parameter file and return a dictionary
    '''

    parameters = {}

    f = file(filename, 'rb')

    for line in f.readlines():
        if '=' in line:
            cols = line.split('=')
            value, key = cols[0].strip(), cols[1].strip()
            parameters[key] = value

    f.close()

    return parameters
