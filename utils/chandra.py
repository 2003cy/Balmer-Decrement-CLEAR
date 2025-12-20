from astropy.table import Table
import pandas as pd
import numpy as np
#functions to load chandra catalogs with custom column specs
def load_chandra_gn(catalog_dir):
    #goodsn xray catalog
    column_specs_gn = [
        (0, 3), (4, 6), (7, 9), (10, 15), (16, 17), (17, 19), (20, 22), (23, 27),
        (28, 33), (34, 36), (37, 40), (41, 46), (47, 55), (56, 62), (63, 69), (70, 78),
        (79, 85), (86, 92), (93, 101), (102, 108), (109, 115), (116, 117), (118, 120),
        (121, 123), (124, 129), (130, 131), (131, 133), (134, 136), (137, 141), (142, 147),
        (148, 153), (154, 161), (162, 164), (165, 167), (168, 173), (174, 175), (175, 177),
        (178, 180), (181, 185), (186, 191), (192, 194), (195, 197), (198, 203), (204, 205),
        (205, 207), (208, 210), (211, 215), (216, 221), (222, 224), (225, 227), (228, 233),
        (234, 235), (235, 237), (238, 240), (241, 245), (246, 251), (252, 254), (255, 257),
        (258, 263), (264, 265), (265, 267), (268, 270), (271, 275), (276, 281), (282, 284),
        (285, 287), (288, 293), (294, 295), (295, 297), (298, 300), (301, 305), (306, 311),
        (312, 314), (315, 317), (318, 323), (324, 325), (325, 327), (328, 330), (331, 335),
        (336, 341), (342, 348), (349, 351), (352, 358), (359, 365), (366, 372), (373, 379),
        (380, 386), (387, 395), (396, 402), (403, 409), (410, 416), (417, 423), (424, 430),
        (431, 437), (438, 444), (445, 447), (448, 450), (451, 456), (457, 458), (458, 460),
        (461, 463), (464, 468), (469, 478), (479, 488), (489, 498), (499, 504), (505, 510),
        (511, 516), (517, 522), (523, 528), (529, 534), (535, 544), (545, 554), (555, 564),
        (565, 574), (575, 581), (582, 585)
    ]
    column_names_gn = [
        "ID", "RAh", "RAm", "RAs", "DE-", "DEd", "DEm", "DEs", "Log(P)", "Log(Pm)",
        "PosErr", "OffAng", "FCts", "UFCts", "LFCts", "SCts", "USCts", "LSCts", "HCts",
        "UHCts", "LHCts", "Ifpsf", "CRAh", "CRAm", "CRAs", "CDE-", "CDEd", "CDEm",
        "CDEs", "C-XOff", "Cmag", "Ccat", "VRAh", "VRAm", "VRAs", "VDE-", "VDEd",
        "VDEm", "VDEs", "Vmag", "GDRAh", "GDRAm", "GDRAs", "GDDE-", "GDDEd", "GDDEm",
        "GDDEs", "GDmag", "CFRAh", "CFRAm", "CFRAs", "CFDE-", "CFDEd", "CFDEm", "CFDEs",
        "CFmag", "KsRAh", "KsRAm", "KsRAs", "KsDE-", "KsDEd", "KsDEm", "KsDEs", "Ksmag",
        "RRAh", "RRAm", "RRAs", "RDE-", "RDEd", "RDEm", "RDEs", "Rmag", "SRAh", "SRAm",
        "SRAs", "SDE-", "SDEd", "SDEm", "SDEs", "Smag", "zspec", "zref", "zphoty",
        "zphotyl", "zphotyu", "q_zphoty", "zphotya", "f_zphoty", "zphots", "zphotsl",
        "zphotsu", "q_zphots", "f_zphots", "zadopt", "A03", "A03RAh", "A03RAm", "A03RAs",
        "A03DE-", "A03DEd", "A03DEm", "A03DEs", "FExp", "SExp", "HExp", "BRat", "UBRat",
        "LBRat", "PInd", "UPInd", "LPInd", "FFlux", "SFlux", "HFlux", "Lx", "Type", "Notes"
    ]

    chandra_gn = Table.from_pandas(
                    pd.read_fwf(f'{catalog_dir}/chandra_gn.txt', 
                                colspecs=column_specs_gn,  
                                names=column_names_gn,     
                                skiprows=217))
    
    #gn has RA DEC in hms dms that we need to convert to degrees
    # Extract RA and DEC columns from the Astropy table
    ra_h = chandra_gn['CFRAh']  # Hour of Right Ascension
    ra_m = chandra_gn['CFRAm']  # Minute of Right Ascension
    ra_s = chandra_gn['CFRAs']  # Second of Right Ascension
    dec_sign = chandra_gn['CFDE-']  # Sign of Declination
    dec_d = chandra_gn['CFDEd']  # Degree of Declination
    dec_m = chandra_gn['CFDEm']  # Arcminute of Declination
    dec_s = chandra_gn['CFDEs']  # Arcsecond of Declination
    # Convert RA to degrees using scipy vectorization
    ra_deg = 15 * (np.add(ra_h, np.add(ra_m / 60, ra_s / 3600)))  # Convert RA to degrees (1 hour = 15 degrees)
    # Convert DEC to degrees based on the sign using scipy functions
    dec_deg = np.add(dec_d, np.add(dec_m / 60, dec_s / 3600))  # Convert Declination to degrees
    dec_deg = np.where(dec_sign == '-', -dec_deg, dec_deg)  # Apply the sign to Declination
    # Add converted RA and DEC in degrees to the table
    chandra_gn['RAdeg'] = ra_deg  # Right Ascension in degrees
    chandra_gn['DEdeg'] = dec_deg  # Declination in degrees
    
    return chandra_gn



def load_chandra_gs(catalog_dir):
    column_specs_gs = [
        (0, 4), (5, 14), (15, 25), (26, 31), (32, 34), (35, 39), (40, 45), 
        (46, 53), (54, 59), (60, 65), (66, 73), (74, 79), (80, 85), (86, 93), 
        (94, 99), (100, 105), (106, 107), (108, 115), (116, 121), (122, 145), 
        (146, 155), (156, 166), (167, 172), (173, 182), (183, 193), (194, 199), 
        (200, 209), (210, 220), (221, 226), (227, 236), (237, 247), (248, 253), 
        (254, 263), (264, 274), (275, 280), (281, 290), (291, 301), (302, 307), 
        (308, 317), (318, 328), (329, 334), (335, 341), (342, 350), (351, 353), 
        (354, 359), (360, 365), (366, 371), (372, 377), (378, 383), (384, 389), 
        (390, 396), (397, 402), (403, 408), (409, 414), (415, 423), (424, 432), 
        (433, 441), (442, 448), (449, 455), (456, 462), (463, 468), (469, 474), 
        (475, 480), (481, 490), (491, 500), (501, 510), (511, 520), (521, 530), 
        (531, 540), (541, 547), (548, 551), (552, 555), (556, 560)
    ]
    column_names_gs = [
        "XID", "RAdeg", "DEdeg", "logPB", "WAV", "PE", "Angle", "FB", "e_FB", "E_FB",
        "SB", "e_SB", "E_SB", "HB", "e_HB", "E_HB", "PNOTE", "CPCAT", "CPDIS", "CPNOTE",
        "WFI-RA", "WFI-DE", "WFImag", "GOODSS-RA", "GOODSS-DE", "GOODSSmag", 
        "GEMS-RA", "GEMS-DE", "GEMSmag", "CANDELS-RA", "CANDELS-DE", "CANDELSmag",
        "TENIS-RA", "TENIS-DE", "TENISmag", "SEDS-RA", "SEDS-DE", "SEDSmag",
        "VLA-RA", "VLA-DE", "VLAmag", "zSpec", "q_zSpec", "r_zSpec", "zL10", "zR11", 
        "zH14", "zS14", "zS15", "zS16", "zF", "r_zF", "e_zF", "E_zF", "FExp", "SExp", 
        "HExp", "BR", "e_BR", "E_BR", "Gamma", "e_Gamma", "E_Gamma", "FFB", "FSB",
        "FHB", "LX", "NH", "LXc", "OType", "X11ID", "X16ID", "R13ID"
    ]

    chandra_gs = Table.from_pandas(
                    pd.read_fwf(f'{catalog_dir}/chandra_gs.txt', 
                                colspecs=column_specs_gs,  
                                names=column_names_gs,     
                                skiprows=174))

    return chandra_gs
