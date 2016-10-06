
# There are only a handful of fortran format codes used by this program,
# so these have been manually converted to python format codes

FORMAT_CODES = {
        '20A4'    : ['{:<4s}']*20,
        'I2,A78'  : ['{:>2d}','{:>78s}'],
        '10I8'    : ['{:>8d}']*10,
        'I8'      : ['{:>8d}'],
        '2I8'     : ['{:>8d}']*2,
        '3I8'     : ['{:>8d}']*3,
        '6F12.7'  : ['{:>12.7F}']*6,
        '3E24.16' : ['{:>24.16E}']*3,
        '3E25.17' : ['{:>25.17E}']*3,
        '5E16.8'  : ['{:>16.8E}']*5,
        'I5,5E15.7'  : ['{:>5d}']+['{:>15.7E}']*5
        }

def fortran_format(fortran_format_code, values):
    format_line = FORMAT_CODES[fortran_format_code]
    num_per_line = len(format_line)
    format_codes = [] if len(values)>0 else ['\n']
    newlines = 0
    while True:
        remaining = len(values) - (len(format_codes) - newlines)
        if remaining <= 0: break
        format_codes.extend( format_line[:min(remaining, num_per_line)] )
        format_codes.append('\n')
        newlines += 1
        
    return ''.join(format_codes).format(*values)
        


