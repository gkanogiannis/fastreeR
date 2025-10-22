#
# fastreeR https://github.com/gkanogiannis/fastreeR
#
# Copyright (C) 2021 Anestis Gkanogiannis <anestis@gkanogiannis.com>
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
#

FROM python:3.10-slim

# Install Java
RUN apt-get update && \
    apt-get install -y openjdk-21-jre && \
    apt-get clean && rm -rf /var/lib/apt/lists/*

# Set working directory
WORKDIR /opt/fastreer

# Copy CLI and backend
COPY fastreeR.py /usr/local/bin/fastreeR
COPY inst/java /opt/fastreer/java

# Let CLI know where JARs are
ENV FASTREER_JAR_DIR=/opt/fastreer/java

# Make CLI executable
RUN chmod +x /usr/local/bin/fastreeR

# Default command
ENTRYPOINT ["fastreeR"]
