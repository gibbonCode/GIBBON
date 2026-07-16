FROM intel/oneapi-toolkit:latest AS oneapi
FROM ubuntu:22.04

ARG DEBIAN_FRONTEND=noninteractive

RUN apt-get update && apt-get install -y --no-install-recommends \
	ca-certificates \
	libfftw3-3 \
	unzip \
	&& rm -rf /var/lib/apt/lists/*

COPY --from=oneapi /opt/intel/oneapi /opt/intel/oneapi

ENV ONEAPI_ROOT=/opt/intel/oneapi
ENV MKLROOT=/opt/intel/oneapi/mkl/latest
ENV LD_LIBRARY_PATH=/opt/intel/oneapi/mkl/latest/lib/intel64:/opt/intel/oneapi/compiler/latest/linux/compiler/lib/intel64_lin:${LD_LIBRARY_PATH}

WORKDIR /tmp/febio

# Private build context must provide the FEBio installer zip.
ARG FEBIO_INSTALLER_ZIP=FEBioStudio_3.1_FEBio_4.12_Linux_installer.zip
COPY ${FEBIO_INSTALLER_ZIP} ./

# RUN set -eux; \
RUN unzip -q "./${FEBIO_INSTALLER_ZIP}" -d ./unpacked
RUN INSTALLER_BIN="$(find ./unpacked -maxdepth 3 -type f -name '*.run' | head -n 1 || true)" && \
	echo $INSTALLER_BIN && \
	chmod u+x "$INSTALLER_BIN" && \
	"$INSTALLER_BIN" --mode unattended --prefix /opt/FEBio
RUN test -x /opt/FEBio/bin/febio4
# RUN /opt/FEBio/bin/febio4 --help >/dev/null

WORKDIR /workspace

RUN rm -rf /tmp/febio

# RUN set -eux; \
# 	unzip -q "./${FEBIO_INSTALLER_ZIP}" -d ./unpacked; \
# 	chmod -R +x ./unpacked; \
# 	INSTALLER_BIN="$(find ./unpacked -maxdepth 3 -type f -name '*.run' | head -n 1 || true)"; \
# 	test -n "$INSTALLER_BIN"; \
# 	"$INSTALLER_BIN" --mode unattended --prefix /opt/FEBio --enable-components comp_febio || "$INSTALLER_BIN" --mode unattended --prefix /opt/FEBio; \
# 	test -x /opt/FEBio/bin/febio4; \
# 	/opt/FEBio/bin/febio4 --help >/dev/null; \
# 	rm -rf /tmp/febio

ENV FEBIO_PATH=/opt/FEBio/bin/febio4
ENV PATH=/opt/FEBio/bin:${PATH}

WORKDIR /workspace
