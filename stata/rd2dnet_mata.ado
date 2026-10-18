*! version 0.2.0  rd2dnet
*! Compile the rd2dnet Mata library (lrd2dnet.mata) once per session.
program define rd2dnet_mata
    version 16
    capture mata: rdn_mata_version()
    if _rc {
        quietly findfile lrd2dnet.mata
        quietly do `"`r(fn)'"'
    }
end
