
# There are only a handful of fortran format codes used by this program,
# so these have been manually converted to python format codes

FORMAT_CODES = {
        '20a4'    : ['{:<4.4s}']*20,
        'i2,a78'  : ['{:>2d}','{:>78.78s}'],
        '10i8'    : ['{:>8d}']*10,
        'i8'      : ['{:>8d}'],
        '2i8'     : ['{:>8d}']*2,
        '3i8'     : ['{:>8d}']*3,
        '6f12.7'  : ['{:>12.7f}']*6,
        '3e24.16' : ['{:>24.16E}']*3,
        '3e25.17' : ['{:>25.17E}']*3,
        '5e16.8'  : ['{:>16.8E}']*5,
        'i5,5e15.7'  : ['{:>5d}']+['{:>15.7E}']*5
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
        


