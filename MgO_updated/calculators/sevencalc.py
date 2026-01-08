from sevenn.calculator import SevenNetCalculator

def sevenn():
    calc = SevenNetCalculator(
        model="7net-omni", 
        modal="mpa",        
        enable_cueq=False, 
        enable_flash=False  
    )
    return calc