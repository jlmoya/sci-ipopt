#!/bin/bash
# macos-fix-arpack-mpi.sh — prerequisite for sci-ipopt on macOS arm64 / Scilab 2027.
#
# WHY: IPOPT's bundled sequential MUMPS (libdmumps/libmumps_common) uses FLAT
# ("dynamically looked up") MPI symbols. Scilab's arnoldi module links Homebrew
# libarpack, which links REAL open-mpi — so loading Scilab force-loads OpenMPI into
# the process. MUMPS's flat MPI symbols then bind to real OpenMPI instead of IPOPT's
# sequential libmpiseq stub → "MPI_Comm_f2c called before MPI_INIT" abort at solve time.
#
# FIX: redirect arpack's 4 open-mpi dependencies to libmpiseq (the sequential stub).
# Scilab uses only SERIAL arpack (dnaupd/dneupd — pure LAPACK, never calls MPI), so the
# stubs are never actually invoked; this just stops real OpenMPI from loading. Each dep
# must point at a DISTINCT file (4 copies of libmpiseq) or dyld errors "duplicate linked
# dylib". Verified safe: spec()/eigs still correct, and ipopt + fmincon then solve.
#
# Re-run after `brew upgrade arpack` (reverts the patch) or `brew reinstall ipopt`
# (removes the shims). Idempotent. Revert: cp the .preipopt backup back + codesign.
set -e
ARP=$(otool -L /opt/homebrew/opt/arpack/lib/libarpack.2.dylib >/dev/null 2>&1 && \
      python3 -c "import os;print(os.path.realpath('/opt/homebrew/opt/arpack/lib/libarpack.2.dylib'))")
SHIMDIR=/opt/homebrew/opt/ipopt/lib
MPISEQ="$SHIMDIR/libmpiseq.dylib"
[ -f "$ARP" ] || { echo "arpack not found"; exit 1; }
[ -f "$MPISEQ" ] || { echo "libmpiseq not found — is ipopt installed? (brew install ipopt)"; exit 1; }

# already patched?
if otool -L "$ARP" 2>/dev/null | grep -q "libmpiseq_shim1"; then
  echo "arpack already patched (libmpiseq shims present)."; exit 0
fi

chmod u+w "$ARP"
[ -f "${ARP}.preipopt" ] || cp -p "$ARP" "${ARP}.preipopt"

for n in 1 2 3 4; do
  cp -pf "$MPISEQ" "$SHIMDIR/libmpiseq_shim$n.dylib"
  install_name_tool -id "$SHIMDIR/libmpiseq_shim$n.dylib" "$SHIMDIR/libmpiseq_shim$n.dylib"
  codesign -f -s - "$SHIMDIR/libmpiseq_shim$n.dylib"
done

i=1
for dep in libmpi_usempif08.40 libmpi_usempi_ignore_tkr.40 libmpi_mpifh.40 libmpi.40; do
  install_name_tool -change "/opt/homebrew/opt/open-mpi/lib/$dep.dylib" \
                            "$SHIMDIR/libmpiseq_shim$i.dylib" "$ARP"
  i=$((i+1))
done
codesign -f -s - "$ARP"
echo "Patched: arpack open-mpi deps redirected to libmpiseq shims. ipopt/fmincon will now solve."
