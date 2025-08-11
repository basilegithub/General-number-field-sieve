# This is the file to write logs in the log file

from datetime import datetime

def format_time(now):
    return str(now.year)+"-"+str(now.month)+"-"+str(now.day)+" "+str(now.hour)+":"+str(now.minute)+":"+str(now.second)

def write_log(logs_path, line):
    now = datetime.now()
    
    log_file = open(logs_path, 'a')
    
    print(format_time(now) + " " + line)
    log_file.write(format_time(now) + " " + line+"\n")
    
    log_file.close()
    
def skip_line(logs_path):
    log_file = open(logs_path, 'a')

    print("")
    log_file.write("\n")
    
    log_file.close()