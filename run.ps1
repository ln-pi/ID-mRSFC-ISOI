# Forward all arguments and the exit code to the project-root Python launcher.
# Activate the desired Conda environment first, or set ISOI_PYTHON explicitly.
$pythonExecutable = if ($env:ISOI_PYTHON) { $env:ISOI_PYTHON } else { 'python' }
& $pythonExecutable (Join-Path $PSScriptRoot 'run.py') @args
exit $LASTEXITCODE
