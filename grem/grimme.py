# -*- coding: utf-8 -*-

#
#    grem Graph EM.
#    Copyright (c) 2021 Be The Match operated by National Marrow Donor Program. All Rights Reserved.
#
#    This library is free software; you can redistribute it and/or modify it
#    under the terms of the GNU Lesser General Public License as published
#    by the Free Software Foundation; either version 3 of the License, or (at
#    your option) any later version.
#
#    This library is distributed in the hope that it will be useful, but WITHOUT
#    ANY WARRANTY; with out even the implied warranty of MERCHANTABILITY or
#    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU Lesser General Public
#    License for more details.
#
#    You should have received a copy of the GNU Lesser General Public License
#    along with this library;  if not, write to the Free Software Foundation,
#    Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307  USA.
#
#    > http://www.fsf.org/licensing/licenses/lgpl.html
#    > http://www.opensource.org/licenses/lgpl-license.php
#

import os
from .EM.run_em import run_em_def


def em(conf_file=""):
    project_dir_in_file = ""
    if conf_file == "":
        use_default_path = True
        conf_file = (
            os.path.dirname(os.path.realpath(__file__))
            + "/conf/minimal-em-configuration.json"
        )
        project_dir_in_file = (
            os.path.dirname(os.path.realpath(__file__)).replace("/grim", "") + "/"
        )

    run_em_def(conf_file, project_dir_in_file=project_dir_in_file)
