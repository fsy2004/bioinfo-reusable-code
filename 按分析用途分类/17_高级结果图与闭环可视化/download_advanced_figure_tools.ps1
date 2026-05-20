$ErrorActionPreference = "Stop"

$root = Split-Path -Parent $MyInvocation.MyCommand.Path
$toolDir = Join-Path $root "external_tools"
New-Item -ItemType Directory -Force -Path $toolDir | Out-Null

$repos = @(
  @{ name = "shapviz"; url = "https://github.com/ModelOriented/shapviz.git" },
  @{ name = "DALEX"; url = "https://github.com/ModelOriented/DALEX.git" },
  @{ name = "shap"; url = "https://github.com/shap/shap.git" },
  @{ name = "cellrank"; url = "https://github.com/theislab/cellrank.git" },
  @{ name = "scvelo"; url = "https://github.com/theislab/scvelo.git" },
  @{ name = "miloR"; url = "https://github.com/MarioniLab/miloR.git" },
  @{ name = "nichenetr"; url = "https://github.com/saeyslab/nichenetr.git" },
  @{ name = "liana"; url = "https://github.com/saezlab/liana.git" },
  @{ name = "pySCENIC"; url = "https://github.com/aertslab/pySCENIC.git" },
  @{ name = "decoupler-py"; url = "https://github.com/scverse/decoupler-py.git" },
  @{ name = "squidpy"; url = "https://github.com/scverse/squidpy.git" },
  @{ name = "cell2location"; url = "https://github.com/BayraktarLab/cell2location.git" },
  @{ name = "Tangram"; url = "https://github.com/broadinstitute/Tangram.git" },
  @{ name = "COMMOT"; url = "https://github.com/zcang/COMMOT.git" },
  @{ name = "pyCirclize"; url = "https://github.com/moshi4/pyCirclize.git" }
)

foreach ($repo in $repos) {
  $dest = Join-Path $toolDir $repo.name
  if (Test-Path $dest) {
    Write-Host "Updating $($repo.name)"
    git -C $dest pull --ff-only
  } else {
    Write-Host "Cloning $($repo.name)"
    git clone --depth 1 $repo.url $dest
  }
}

Write-Host "Done. Tools are under $toolDir"

