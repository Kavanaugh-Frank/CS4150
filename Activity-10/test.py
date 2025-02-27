features = [
    "Hist1", "Vmn", "LAD", 
    "RNAPII-S2P", "RNAPII-S5P", "RNAPII-S7P", 
    "Enhancer", "H3K9me3", "H3K20me3", "h3k27me3", "H3K36me3", 
    "NANOG", "pou5f1", "sox2", "CTCF-7BWU"
]

feature_map = {}
with open("Hist1_region_features.csv", "r") as f:
    headers = f.readline().split(",")
    headers = [header.strip() for header in headers]
    print(headers)
    for feature in features:
        feature_map[feature] = headers.index(feature)
print(feature_map)
    
