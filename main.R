
root <- rprojroot::find_root(rprojroot::has_file(".gitignore"))

# create folders if needed
if (dir.exists(file.path(root, "tables"))) dir.create(file.path(root, "tables"))
if (dir.exists(file.path(root, "figures"))) dir.create(file.path(root, "figures"))
if (dir.exists(file.path(root, "data"))) dir.create(file.path(root, "data"))

# Cohort generation
source(file.path(root, "scripts", "summary-stats", "cohort-gen.R"))

# Cohort generation
source(file.path(root, "scripts", "summary-stats", "patient-tables.R"))

# Process of Care
source(file.path(root, "scripts", "summary-stats", "process-of-care.R"))

# Missingness Table
source(file.path(root, "scripts", "summary-stats", "missingness.R"))

# Univariate Analyses
source(file.path(root, "scripts", "univariate", "univariate-bmi.R"))
source(file.path(root, "scripts", "univariate", "univariate-cond.R"))
source(file.path(root, "scripts", "univariate", "association.R"))

# Multivariate Analyses
source(file.path(root, "scripts", "multivariate", "logistic.R"))
