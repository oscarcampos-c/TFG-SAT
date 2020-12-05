param($jsonfile,$resultfilename,$blobfilename)
$json = Get-Content -Raw -Path $jsonfile | ConvertFrom-Json
$connectionstring = ""
$json | foreach {
    $name = $_.name
    $_.properties | foreach {
        if($_.type -eq "AzureBlobStorage")
		{
			$connectionstring = $_.typeProperties.connectionString
		}
    }
}
$containerName = "data"
$ctx = New-AzStorageContext -ConnectionString $connectionstring
Set-AzStorageBlobContent -File "$resultfilename"  -Container $containerName  -Blob "$blobfilename"  -Context $ctx -Force

