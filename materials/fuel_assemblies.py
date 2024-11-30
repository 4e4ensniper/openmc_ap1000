fa_types = [
    {"name": "Z40",
     "enrichment": 4.0,
     "grey_enrichment": 0,
     "gdo2_wo": 0,
     "grey_pos": []
    },
    {"name": "Z13",
     "enrichment": 1.3,
     "grey_enrichment": 0,
     "gdo2_wo": 0,
     "grey_pos": []
    },
    {"name": "Z24",
     "enrichment": 2.4,
     "grey_enrichment": 0,
     "gdo2_wo": 0,
     "grey_pos": []
    },
    {"name": "Z40D2",
     "enrichment": 4.0,
     "grey_enrichment": 3.3,
     "gdo2_wo": 5.0,
     "grey_pos": [114, 118, 122, 126, 130, 134, 138, 142, 146, 150, 154, 158]
    },
    {"name": "Z40D6",
     "enrichment": 4.0,
     "grey_enrichment": 3.3,
     "gdo2_wo": 5.0,
     "grey_pos": [114, 122, 130, 138, 146, 154]
    },
    {"name": "Z40E9",
     "enrichment": 4.0,
     "grey_enrichment": 3.3,
     "gdo2_wo": 5.0,
     "grey_pos": [114, 122, 130, 138, 146, 154, 271, 279, 287]
    },
    {"name": "Z40U1",
     "enrichment": 4.0,
     "grey_enrichment": 3.3,
     "gdo2_wo": 8.0,
     "grey_pos": [114, 122, 130, 138, 146, 154,
                  164, 167, 171, 174, 178, 181, 185, 188, 192, 195, 199, 202,
                  296, 302, 308]
    },
    {"name": "Z40U5",
     "enrichment": 4.0,
     "grey_enrichment": 3.3,
     "gdo2_wo": 8.0,
     "grey_pos": [114, 118, 122, 126, 130, 134, 138, 142, 146, 150, 154, 158,
                  271, 279, 287]
    },
    {"name": "Z44A6",
     "enrichment": 4.4,
     "grey_enrichment": 3.6,
     "gdo2_wo": 5.0,
     "grey_pos": [271, 275, 279, 283, 287, 291]
    },
    {"name": "Z44B2",
     "enrichment": 4.4,
     "grey_enrichment": 3.6,
     "gdo2_wo": 5.0,
     "grey_pos": [114, 118, 122, 126, 130, 134, 138, 142, 146, 150, 154, 158]
    },
    {"name": "Z44Y7",
     "enrichment": 4.4,
     "grey_enrichment": 3.6,
     "gdo2_wo": 8.0,
     "grey_pos": [63, 66, 72, 75, 81, 84, 90, 93, 99, 102, 108, 111,
                   114, 122, 130, 138, 146, 154,
                   205, 211, 217, 223, 229, 235,
                   312, 316, 320]
    },
    {"name": "Z49",
     "enrichment": 4.95,
     "grey_enrichment": 0,
     "gdo2_wo": 0,
     "grey_pos": []
    },
    {"name": "Z49A2",
     "enrichment": 4.95,
     "grey_enrichment": 3.6,
     "gdo2_wo": 5.0,
     "grey_pos": [114, 122, 130, 138, 146, 154, 271, 275, 279, 283, 287, 291]
    },
    {"name": "Z33Z2",
     "enrichment": 3.3,
     "grey_enrichment": 2.4,
     "gdo2_wo": 8.0,
     "grey_pos": [114, 122, 130, 138, 146, 154, 271, 275, 279, 283, 287, 291]
    },
    {"name": "Z49B6",
     "enrichment": 4.95,
     "grey_enrichment": 3.6,
     "gdo2_wo": 5.0,
     "grey_pos": [114, 122, 130, 138, 146, 154]
    },
    {"name": "Z49B9",
     "enrichment": 4.95,
     "grey_enrichment": 3.6,
     "gdo2_wo": 5.0,
     "grey_pos": [114, 118, 122, 130, 134, 138, 146, 150, 154]
    },
    {"name": "Z33Z9",
     "enrichment": 3.3,
     "grey_enrichment": 2.4,
     "gdo2_wo": 8,
     "grey_pos": [114, 118, 122, 130, 134, 138, 146, 150, 154]############
    },
    {"name": "Z49Y7",
     "enrichment": 4.95,
     "grey_enrichment": 3.6,
     "gdo2_wo": 8.0,
     "grey_pos": [63, 66, 72, 75, 81, 84, 90, 93, 99, 102, 108, 111,
                   114, 122, 130, 138, 146, 154,
                   205, 211, 217, 223, 229, 235,
                   312, 316, 320]
    },
    {"name": "Z49Y9",
     "enrichment": 4.95,
     "grey_enrichment": 3.6,
     "gdo2_wo": 8.0,
     "grey_pos": [114, 122, 130, 138, 146, 154, 271, 279, 287]
    }
]

def find_name(name, fas):
    for fa in fas:
        if fa["name"] == name:
            return fa
    return None
