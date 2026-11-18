param([string]$File)
$here = Split-Path -Parent $MyInvocation.MyCommand.Path
Add-Type -Path "$here\Msaa.cs" -ReferencedAssemblies System.Drawing, Accessibility, System.Windows.Forms
$p = Get-Process "RASAero II" -ErrorAction Stop | Select-Object -First 1
$dlg = $null; $main = $null
foreach ($w in [Msaa]::TopWindows($p.Id)) { if ([Msaa]::Class($w) -eq '#32770') { $dlg = $w } else { $main = $w } }
if (-not $dlg) {
  $known = New-Object 'System.Collections.Generic.List[IntPtr]'
  foreach ($w in [Msaa]::TopWindows($p.Id)) { $known.Add($w) }
  [void][Drive]::PressAsync($main, "Open")
  $dlg = [Drive]::WaitNew($p.Id, $known, 8000)
}
$edit = [Drive]::Child($dlg, 1148)
if ([Msaa]::Class($edit) -eq 'ComboBoxEx32') { foreach ($c in [Msaa]::Children($edit)) { if ([Msaa]::Class($c) -eq 'Edit') { $edit = $c } } }
[Msaa]::SetText($edit, $File)
Start-Sleep -Milliseconds 300
[Msaa]::Click([Drive]::Child($dlg, 1))
Start-Sleep -Seconds 3
foreach ($w in [Msaa]::TopWindows($p.Id)) { "window: '" + [Msaa]::Text($w) + "' class " + [Msaa]::Class($w) }
[Msaa]::DumpAcc($main, 6, $false) -split "`n" | Where-Object { $_ -match "All Turbulent" }
