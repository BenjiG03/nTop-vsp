param([string]$OutDir, [string]$Tag = "", [string]$Alphas = "0,2,4,8")
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
Add-Type -Path "$here\Msaa.cs" -ReferencedAssemblies System.Drawing, Accessibility, System.Windows.Forms
$p = Get-Process "RASAero II" -ErrorAction Stop | Select-Object -First 1

function Windows { , @([Msaa]::TopWindows($p.Id)) }
function Find-RunTest {
  foreach ($w in [Msaa]::TopWindows($p.Id)) { if ([Msaa]::Text($w) -eq 'Run Test') { return $w } }
  return [IntPtr]::Zero
}
function Main-Window {
  foreach ($w in [Msaa]::TopWindows($p.Id)) { if ([Msaa]::Text($w) -like 'RASAero II*') { return $w } }
  throw "no main window"
}

foreach ($alpha in ($Alphas -split ",")) {
  $dlg = Find-RunTest
  if ($dlg -eq [IntPtr]::Zero) {
    $known = New-Object 'System.Collections.Generic.List[IntPtr]'
    foreach ($w in [Msaa]::TopWindows($p.Id)) { $known.Add($w) }
    [void][Drive]::PressAsync((Main-Window), "Run Test")
    $dlg = [Drive]::WaitNew($p.Id, $known, 8000)
    Start-Sleep -Milliseconds 500
  }
  $fileEdit = $null; $alphaEdit = $null; $button = $null
  foreach ($c in [Msaa]::Children($dlg)) {
    $cls = [Msaa]::Class($c); $t = [Msaa]::Text($c)
    if ($cls -like '*.EDIT.*') {
      if ($t -match '[:\\]') { $fileEdit = $c }
      elseif ($t -notmatch '^0\.0000$') { $alphaEdit = $c }
    }
    if ($cls -like '*.BUTTON.*' -and $t -eq 'Run Test') { $button = $c }
  }
  if (-not ($fileEdit -and $alphaEdit -and $button)) { throw "controls not found: $([Msaa]::Describe($dlg))" }
  $out = Join-Path $OutDir ("runtest{0}_alpha{1}.txt" -f $Tag, $alpha)
  if (Test-Path $out) { Remove-Item $out -Confirm:$false }
  [Msaa]::SetText($alphaEdit, ("{0}" -f $alpha))
  [Msaa]::SetText($fileEdit, $out)
  Start-Sleep -Milliseconds 300
  $known = New-Object 'System.Collections.Generic.List[IntPtr]'
  foreach ($w in [Msaa]::TopWindows($p.Id)) { $known.Add($w) }
  [Msaa]::Click($button)
  $sw = [Diagnostics.Stopwatch]::StartNew()
  $done = $false; $last = -1; $stable = 0
  while ($sw.Elapsed.TotalSeconds -lt 120) {
    Start-Sleep -Milliseconds 300
    # A message box or a replace prompt would stop the run; report it and accept it.
    foreach ($w in [Msaa]::TopWindows($p.Id)) {
      if (-not $known.Contains($w) -and [Msaa]::Class($w) -eq '#32770') {
        "  popup: '$([Msaa]::Text($w))' :: " + (([Msaa]::Children($w) | ForEach-Object { [Msaa]::Text($_) }) -join ' | ')
        $ok = [Drive]::Child($w, 1)
        if ($ok -ne [IntPtr]::Zero) { [Msaa]::Click($ok) } else { [Msaa]::Close($w) }
        Start-Sleep -Milliseconds 300
      }
    }
    if (Test-Path $out) {
      $now = (Get-Item $out).Length
      if ($now -gt 100000 -and $now -eq $last) { $stable++ } else { $stable = 0 }
      $last = $now
      if ($stable -ge 6) { $done = $true; break }
    }
  }
  $size = if (Test-Path $out) { (Get-Item $out).Length } else { -1 }
  "alpha $alpha done=$done size=$size t=$([int]$sw.Elapsed.TotalSeconds)s dialog-open=$((Find-RunTest) -ne [IntPtr]::Zero)"
  $still = Find-RunTest
  if ($still -ne [IntPtr]::Zero -and $size -gt 0) { [Msaa]::Close($still); Start-Sleep -Milliseconds 500 }
}
