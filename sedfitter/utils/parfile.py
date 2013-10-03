from __future__ import print_function, division


def read(filename, format):
    """
    Read a parameter file and return a dictionary
    """

    parameters = {}

    for line in open(filename):
        if '=' in line and line[0] != "#" and line.strip() != "":
            cols = line.split('=')
            if format == 'par':
                value, key = cols[0].strip(), cols[1].strip()
            elif format == 'conf':
                key, value = cols[0].strip(), cols[1].strip()
            else:
                raise Exception("Format should be par or conf")
            try:
                value = int(value)
            except:
                try:
                    value = float(value)
                except:
                    if value.lower() in ['y', 'yes']:
                        value = True
                    elif value.lower() in ['n', 'no']:
                        value = False
            parameters[key] = value

    return parameters
