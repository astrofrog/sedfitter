def read(filename):
    
    f = file(filename,'rb')
    
    parameters = {}
    
    for line in f.readlines():
        if '=' in line:
            cols = line.split('=')
            value,key = cols[0].strip(),cols[1].strip()
            parameters[key] = value
            
    f.close()    
        
    return parameters